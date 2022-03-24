open Optitrust
open Target
open Ast

let add_prefix (prefix : string) (indices : string list) : string list =
    List.map (fun x -> prefix ^ x) indices

let step = cTopFunDef "step"
let stepLF = cTopFunDef "stepLeapFrog"
let stepsl = [stepLF; step]
let repPart = cTopFunDef "reportParticlesState"
let addPart = cTopFunDef "addParticle"

let dims = ["X"; "Y"; "Z"]
let nb_dims = List.length dims
let iter_dims f = List.iter f dims
let map_dims f = List.map f dims
let idims = map_dims (fun d -> "i" ^ d)
let delocalize_sum = Local_arith (Lit_double 0., Binop_add)
let delocalize_bag = Local_obj ("bag_init_initial", "bag_append", "bag_free_initial")
let align = 64

(* Grab the "usechecker" flag from the command line *)
let usechecker = ref false
let _= Run.process_cmdline_args
  [("-usechecker", Arg.Set usechecker, " use -DCHECKER as preprocessor flag")]
  (* LATER: use a generic -D flag for optitrust *)
let usechecker = !usechecker


(* UNCOMMENT THIS LINE FOR WORKING ON THE VERSION WITH THE CHECKER
let usechecker = true
*)

let onlychecker p = if usechecker then [p] else []
let doublepos = false (* LATER: Arthur make this a command line command *)
let doublepos = if usechecker then false else doublepos

let stepFuns =
  (if usechecker then [repPart] else [])
     @ stepsl

let stepsReal = cOr (List.map (fun f -> [f]) stepsl) (* LATER: rename *)
let steps = cOr (List.map (fun f -> [f]) stepFuns)

let prepro = onlychecker "-DCHECKER"
let prepro = ["-DPRINTPERF"] @ prepro

(* LATER let prefix = if usechecker then "pic_demo_checker" else "pic_demo"
   ~prefix *)

let _ = Run.script_cpp ~parser:Parsers.Menhir ~prepro ~inline:["pic_demo.h";"bag.hc";"particle.hc";"optitrust.h";"bag_atomics.h";"bag.h-"] (fun () ->

  Printf.printf "CHECKER=%d\n" (if usechecker then 1 else 0);

  (* Part 1: sequential optimizations *)

  bigstep "Optimization and inlining of [matrix_vect_mul]";
  let ctx = cTopFunDef "matrix_vect_mul" in
  !! Function.inline [ctx; cOr [[cFun "vect_mul"]; [cFun "vect_add"]]];
  !! Struct.set_explicit [nbMulti; ctx; cWriteVar "res"];
  !! Loop.fission ~split_between:true [ctx; cFor "k"];
  !! Loop.unroll [nbMulti; ctx; cFor "k"];
  !! Instr.accumulate ~nb:8 [nbMulti; ctx; sInstrRegexp ~substr:true "res.*\\[0\\]"];
  !! Function.inline ~delete:true [nbMulti;cFun "matrix_vect_mul"];

  bigstep "Optimization in [cornerInterpolationCoeff], before it is inlined";
  let ctx = cTopFunDef "cornerInterpolationCoeff" in
  let ctx_rv = cChain [ctx; sInstr "r.v"] in (* LATER rewrite pour cX *)
  (* BEAUTIFY
    in pic_demo.c
    in double_nbCorners cornerInterpolationCoeff(vect pos)
      const double cX = 1. + -1. * rX;
    should be replaced with
      const double cX = 1. - rX;
    and we should use a rewrite rule to obtain the previous code.
  *)
  !! Rewrite.equiv_at "double a; ==> a == (0. + 1. * a)" [nbMulti; ctx_rv; cVar ~regexp:true "r."];
  !! Variable.inline [nbMulti; ctx; cVarDef ~regexp:true "c."];
  !! Variable.intro_pattern_array ~const:true ~pattern_aux_vars:"double rX, rY, rZ"
      ~pattern_vars:"double coefX, signX, coefY, signY, coefZ, signZ"
      ~pattern:"(coefX + signX * rX) * (coefY + signY * rY) * (coefZ + signZ * rZ)"
      [nbMulti; ctx_rv; dRHS];
  !! Instr.move_out ~dest:[tBefore; ctx] [nbMulti; ctx; cVarDef ~regexp:true "\\(coef\\|sign\\)."];
  !! Loop.fold_instrs ~index:"k" [cTopFunDef "cornerInterpolationCoeff"; cCellWrite ~base:[cVar "r"] ()];

  bigstep "Update particles in-place instead of in a local variable ";
  !! Variable.reuse ~space:(expr "p->speed") [step; cVarDef "speed2" ];
  !! Variable.reuse ~space:(expr "p->pos") [step; cVarDef "pos2"];

  bigstep "Reveal write operations involved in the manipulation of particles and vectors";
  let ctx = cOr [[cFunDef "bag_push_serial"]; [cFunDef "bag_push_concurrent"]] in
  !! List.iter (fun typ -> Struct.set_explicit [nbMulti; ctx; cWrite ~typ ()]) ["particle"; "vect"];
  !! Function.inline [steps; cOr [[cFun "vect_mul"]; [cFun "vect_add"]]];
  !! Trace.reparse ();
  !! Struct.set_explicit [nbMulti; steps; cFieldWrite ~base:[cVar "p"] ()];

  bigstep "Inlining of [cornerInterpolationCoeff] and [accumulateChargeAtCorners]";
  !! Function.inline [nbMulti; cFunDef "cornerInterpolationCoeff"; cFun ~regexp:true "relativePos."];
  !! Function.inline [step; cFun "accumulateChargeAtCorners"];
  !! Function.inline ~vars:(AddSuffix "2") [step; cFun "idCellOfPos"];
  !! List.iter (fun f -> Function.inline ~vars:(AddSuffix "${occ}") [nbMulti; f; cFun "cornerInterpolationCoeff"])
     stepsl;
  !! iter_dims (fun d -> Variable.reuse ~space:(var ("i" ^ d ^ "2")) [step; cVarDef ("i" ^ d ^ "1")]);

  bigstep "Optimization of charge accumulation";
  !! Sequence.intro ~mark:"fuse" ~start:[step; cVarDef "contribs"] ();
  !! Loop.fusion_targets [cMark "fuse"];
  !! Trace.reparse();
  !! Instr.inline_last_write [step; cCellRead ~base:[cFieldRead ~base:[cVar "contribs"] ()] ()];

  bigstep "Low level iteration on chunks of particles";
  (* Since function inline transformation was not designed to work with calls inside for loops, here we need
      first to bind a variable to that call then inline that function call and finally inline the binding variable *)
  !! Function.bind_intro ~fresh_name:"bag_iter" [steps; cOr [[cFun "bag_iter_begin"]; [cFun "bag_iter_destructive_begin"]]];
  !! Function.inline [steps; cOr [[cFun "bag_iter_begin"]; [cFun "bag_iter_destructive_begin"]]];
  !! Variable.inline [steps; cVarDef "bag_iter"];

  !! Sequence.intro ~mark:"loop" ~start:[steps; cVarDef "bag_it"] ~nb:2 ();
  !! Sequence.intro_on_instr [steps; cMark "loop"; cFor_c ""; dBody];
  !! Function_basic.uninline ~fct:[cFunDef "bag_iter_ho_basic"~body:[cVarDef "it"]] [steps; cMark "loop"]; (* TODO Calling Function.uninline loops *)
  !! Expr.replace_fun "bag_iter_ho_chunk" [steps; cFun "bag_iter_ho_basic"];
  !! Function.inline [steps; cFun "bag_iter_ho_chunk"];
  !! List.iter (fun f -> Function.beta ~indepth:true [f]) stepFuns;
  !! Instr.delete [nbMulti; cTopFunDef ~regexp:true "bag_iter.*"];

  bigstep "Elimination of the pointer on a particle, to prepare for aos-to-soa";
  !! Variable.init_detach [steps; cVarDef "p"];
  !! Function.inline ~delete:true ~vars:(AddSuffix "${occ}") [nbMulti; step; cFun "wrapArea"]; (* BEAUTIFY: move elsewhere *)
  !! Variable.inline [nbMulti; step; cVarDef ~regexp:true "[xyz]."];  (* BEAUTIFY: move elsewhere? *)
  !! Instr.inline_last_write [nbMulti; steps; cRead ~addr:[cStrictNew; cVar "p"] ()];
  !! Instr.delete [steps; cVarDef "p"];

  bigstep "AOS-TO-SOA";
  (* BEAUTIFY: discuss how the first 2 operations could be gathered with the
     "Reveal write operations involved" step *)
  !! Struct.set_explicit [step; cVarDef "p2"];
  !! Struct.set_explicit [nbMulti; step; cFieldWrite ~base:[cVar "p2"] ~regexp:true ~field:"\\(speed\\|pos\\)" ()];
  !! List.iter (fun f -> Struct.inline f [cTypDef "particle"]) ["speed"; "pos"];
  !! Struct.inline "items" [cTypDef "chunk"];

  bigstep "Apply scaling factors on the electric field";
  !! Struct.to_variables [steps; cVarDef "fieldAtPos"];
  !! Variable.insert_list_same_type ~reparse:true (atyp "const double") (["factorC", expr "particleCharge * stepDuration * stepDuration / particleMass"]
      @ (map_dims (fun d -> ("factor" ^ d, expr ("factorC / cell" ^ d))))) [occFirst; tBefore; step; cFor "idCell"];
  !! Function.inline ~delete:true [steps; cFun "getFieldAtCorners"];
  !! Variable.rename ~into:"field_at_corners" [step; cVarDef "res"];
  !! Struct.set_explicit [step; cFor "k"; cCellWrite ~base:[cFieldRead ~base:[cVar "field_at_corners"] ()] ()];
  !! iter_dims (fun d ->
      let d1 = String.lowercase_ascii d in
      Accesses.scale ~factor:(var ("factor" ^ d)) [step; cFor "k"; cFieldWrite ~field:d1 ()]);
  !! iter_dims (fun d ->
       Accesses.scale ~factor:(var ("factor" ^ d)) [step; cVarDef "accel"; cReadVar ("fieldAtPos" ^ d)]);
  !! Variable.unfold [step; cVarDef  "factorC"];
  !! Variable.unfold ~at:[cVarDef "accel"] [nbMulti; step; cVarDef ~regexp:true "factor."];
  !! Arith.(simpl ~indepth:true expand) [nbMulti; step; cVarDef "accel"];

  bigstep "Applying a scaling factor on speeds";
  !! Struct.set_explicit [addPart; cVarDef "p"];
  !! iter_dims (fun d ->
      let d_lc = String.lowercase_ascii d in
      Accesses.scale ~factor:(expr ("(cell"^d^"/stepDuration)")) [addPart; cFieldRead ~field:d_lc ~base:[cVar "speed"] ()]);
  !! iter_dims (fun d ->
      Accesses.scale ~factor:(expr ("(stepDuration / cell"^d^")"))
      [nbMulti; steps; cFieldWrite ~base:[cVar "c"] (); sExprRegexp ~substr:true ("c->itemsSpeed" ^ d ^ "\\[i\\]")]);
  if usechecker then (!! iter_dims (fun d ->
        Accesses.scale ~neg:true ~factor:(expr ("(cell"^d^"/stepDuration)")) [repPart; cVarDef ("speed"^d); cInit()]));

  bigstep "Applying a scaling factor on positions";
  !! iter_dims (fun d ->
      let d_lc = String.lowercase_ascii d in
      Accesses.scale ~factor:(expr ("cell"^d)) [addPart; cFieldRead ~field:d_lc ~base:[cVar "pos"] ()]);
  !! iter_dims (fun d ->
     Accesses.scale ~neg:true ~factor:(expr ("cell"^d))
         [nbMulti; steps; cOr [[sExprRegexp ~substr:true ("c->itemsPos" ^ d ^ "\\[i\\]")]; [cFieldWrite ~field:("pos"^d)()]]]); (* BEAUTIFY: why is the subst:true needed above? *)

  bigstep "Simplify arithmetic expressions after scaling";
  !! Trace.reparse();
  !! Variable.inline [steps; cVarDef "accel"];
  !! Arith.with_nosimpl [nbMulti; steps; cFor "k"] (fun () ->
       Arith.(simpl ~indepth:true expand) [nbMulti; steps]);
       (* BEAUTIFY: remind me why we can't do a infixop just here? it would be very pretty;
          (even if we have to undo it later) *)

  bigstep "Enumerate grid cells by coordinates";
  (* BEAUTIFY: inline nbCells, then call grid_enumerate, which should parse the product automatically i<X*Y*Z *)
  !! Loop.grid_enumerate (map_dims (fun d -> ("i" ^ d, "grid" ^ d))) [step; cFor "idCell" ~body:[cFor "k"]];

  bigstep "Code cleanup in preparation for shifting of positions";
  !! iter_dims (fun d ->
      Variable.bind ~const:true ~typ:(Some (atyp "double")) ("p" ^ d ^ "2") [occLast;step; cCellWrite ~base:[cFieldRead ~field:("itemsPos" ^ d) ()] (); dRHS];
      Variable.bind ~const:true ("p" ^ d ) ~typ:(Some (atyp "double")) [occFirst;step; cCellWrite ~base:[cFieldRead ~field:("itemsPos" ^ d) ()] (); dRHS]);
  !! iter_dims (fun d ->
      Instr.inline_last_write [step; cVarDef ~regexp:true "p.2"; cCellRead ~base:[cFieldRead ~field:("itemsPos" ^ d) ()]()]);
  !! iter_dims (fun d ->
      Instr.read_last_write [step; cVarDef ~regexp:true "i.2"; cCellRead ~base:[cFieldRead ~field:("itemsPos" ^ d) ()]()]);
  !! Instr.(gather_targets ~dest:GatherAtFirst) [step; cVarDef ~regexp:true ("\\(i.*2\\|p.2\\|i.2\\)")];
  !! Instr.(gather_targets ~dest:(GatherAt [tBefore; cVarDef "p2"])) [step; cVarDef ~regexp:true "r.1"];

  bigstep "Shifting of positions: make positions relative to the containing cell";
  !! Instr.move ~dest:[tBefore; addPart; cVarDef "p"] [addPart; cVarDef "idCell"];
  !! List.iter (fun tg ->
      Variable.insert ~typ:(atyp "coord") ~name:"co" ~value:(expr "coordOfCell(idCell)") tg)
      ([ [tAfter; addPart; cVarDef "idCell"] ] @ onlychecker [tFirst; repPart; cFor "idCell"; dBody]); (* LATER: make cOr work for targetBetweens (hard) *)
  !! iter_dims (fun d ->
      Accesses.shift ~neg:true ~factor:(expr ("co.i"^d)) [cOr (
        [[addPart; cFieldWrite ~field:("pos"^d) ()]] @
        (onlychecker [repPart; cCellRead ~base:[cFieldRead ~field:("itemsPos" ^ d) ()] ()]) )] );
  !! iter_dims (fun d ->
      Accesses.shift ~neg:true ~factor:(expr ("i" ^ d ^ "0")) [stepsReal; cVarDef ~regexp:true "r.0"; cCellRead ~base:[cFieldRead ~field:("itemsPos" ^ d) ()]()];
      Accesses.shift ~neg:true ~factor:(expr ("i" ^ d)) [step; cVarDef ~regexp:true ("p" ^ d); cCellRead ~base:[cFieldRead ~field:("itemsPos" ^ d) ()]()];
      Accesses.shift ~neg:true ~factor:(expr ("i" ^ d ^ "2")) [step; cCellWrite ~base:[cFieldRead ~field:("itemsPos" ^ d) ()]()];
      Accesses.shift ~neg:true ~factor:(expr ("i" ^ d ^ "2")) [step; cVarDef ("r" ^ d ^ "1"); cCellRead ~base:[cFieldRead ~field:("itemsPos" ^ d) ()]()]);

  bigstep "Simplify arithmetic expressions after shifting of positions";
  !! Trace.reparse();
  !! Arith.with_nosimpl [nbMulti; stepsReal; cFor "k"] (fun () ->
       Arith.(simpl ~indepth:true expand) [nbMulti; stepsReal]);
  !! Instr.delete [nbMulti; stepsReal; cVarDef ~regexp:true "i.0"];
  !! Variable.fold ~at:[cFieldWrite ~base:[cVar "p2"] ()] [nbMulti; step; cVarDef ~regexp:true "r.1"];
  !! Rewrite.equiv_at ~glob_defs:"double fwrap(double, double);\n" "double x, y, z; ==> (fwrap(x,y)/z) == (fwrap(x/z, y/z))" [cVarDef ~regexp:true "p.2"; cInit()];
  !! iter_dims (fun d ->
      Instr.read_last_write ~write:[cTopFunDef "computeConstants"; cWriteVar ("cell"^d)] [nbMulti;step; cFun "fwrap";cReadVar ("cell"^d)];);
  !! Arith.(simpl_rec expand) [nbMulti; step; cVarDef ~regexp:true "p.2"];
  (* BEAUTIFY: can we factorize the two calls to simpl_rec in this big step ?*)

  if doublepos then begin
    bigstep "Turn positions into floats, decreasing precision but allowing to fit a larger number of particles";
    !! Cast.insert (atyp "float") [sExprRegexp ~substr:true "p.2 - i.2"];
    !! Struct.update_fields_type "itemsPos." (atyp "float") [cTypDef "chunk"];
  end;

  bigstep "Replacement of the wrap-around operation on doubles with a bitwise operation, assuming grid sizes to be powers of 2";
   let fwrapInt = "double fwrapInt(int m, double v) {
      const int q = int_of_double(v);
      const double r = v - q;
      const int j = wrap(m, q);
      return j + r;
    }" in
  !! Sequence.insert ~reparse:true (stmt fwrapInt) [tBefore; step];
  !! Expr.replace_fun "fwrapInt" [nbMulti;step; cFun "fwrap"];
  !! iter_dims (fun d ->
      Function.inline ~vars:(AddSuffix d) [step; cVarDef ("p"^d^"2"); cFun "fwrapInt"]);
      (* BEAUTIFY the next two steps need to be performed differently, see the paper *)
  !! Instr.delete [nbMulti; step; cVarDef ~regexp:true "i.2"];
  !! iter_dims (fun d ->
      Variable.rename ~into:("i"^d^"2") [step;cVarDef ("j"^d)];);
  !! Variable.inline [nbMulti; step; cVarDef ~regexp:true "p.2"];
  !! Arith.(simpl_rec expand) [nbMulti; step; cCellWrite ~base:[cFieldRead ~regexp:true ~field:("itemsPos.") ()] ()];
  !! Expr.replace_fun "wrapPowerof2" [nbMulti; step; cFun "wrap"];
  !! Function.inline [nbMulti; step; cFun "wrapPowerof2"];

  bigstep "Introduce matrix operations, and prepare loop on charge deposit";
  !! Label.add "core" [step; cFor "iX" ];
  !! Matrix_basic.intro_mmalloc [nbMulti; cFunDef "allocateStructures";cFun "malloc"];
  !! Matrix.intro_mindex (expr "nbCells") [nbMulti; step; cCellAccess ~base:[cOr [[cVar "deposit"]; [cVar "bagsNext"]]]() ];
  !! Label.add "charge" [step; cFor "k" ~body:[cVar "deposit"]];
  !! Variable.inline [occLast; step; cVarDef "indices"];

  bigstep "Duplicate the charge of a corner for the 8 surrounding cells";
  let alloc_instr = [cFunDef "allocateStructures"; cWriteVar "deposit"] in
  !! Matrix.delocalize "deposit" ~into:"depositCorners" ~last:true ~indices:["idCell"] ~init_zero:true
     (* ~labels:["alloc"; ""; "dealloc"] ~dealloc_tg:(Some [cTopFunDef ~regexp:true "dealloc.*"; cFor ""]) *)
     ~dim:(expr "8") ~index:"idCorner" ~acc:"sum" ~ops:delocalize_sum ~use:(Some (expr "k")) ~alloc_instr [cLabel "core"];

  bigstep "Apply a bijection on the array storing charge to vectorize charge deposit";
  let mybij_def =
      "int mybij(int nbCells, int nbCorners, int idCell, int idCorner) {
        coord coord = coordOfCell(idCell);
        int iX = coord.iX;
        int iY = coord.iY;
        int iZ = coord.iZ;
        int res[8] = {
          cellOfCoord(iX, iY, iZ),
          cellOfCoord(iX, iY, wrap(gridZ,iZ-1)),
          cellOfCoord(iX, wrap(gridY,iY-1), iZ),
          cellOfCoord(iX, wrap(gridY,iY-1), wrap(gridZ,iZ-1)),
          cellOfCoord(wrap(gridX,iX-1), iY, iZ),
          cellOfCoord(wrap(gridX,iX-1), iY, wrap(gridZ,iZ-1)),
          cellOfCoord(wrap(gridX,iX-1), wrap(gridY,iY-1), iZ),
          cellOfCoord(wrap(gridX,iX-1), wrap(gridY,iY-1), wrap(gridZ,iZ-1)),
        };
      return MINDEX2(nbCells, nbCorners, res[idCorner], idCorner);
      }" in
  !! Sequence.insert (stmt mybij_def) [tBefore; step];
  !! Matrix.biject "mybij" [step; cVarDef "depositCorners"];
  !! Expr.replace ~reparse:false (expr "MINDEX2(nbCells, 8, idCell2, k)")
      [step; cLabel "charge"; cFun "mybij"];

  (* Part 2: parallelization *)

  bigstep "Decompose the loop to allow for parallelization per blocks";
  !! Variable.insert_list_same_type (atyp "const int") [("block", lit "2"); ("halfBlock", (lit "1"))] [tBefore; cVarDef "nbCells"];
  let colorize (tile : string) (color : string) (d:string) : unit =
    let bd = "b" ^ d in
    Loop.tile tile ~bound:TileBoundDivides ~index:("b"^d) [step; cFor ("i"^d)];
    Loop.color (expr color) ~index:("c"^d) [step; cFor bd]
    in
  !! iter_dims (fun d -> colorize "2" "block" d);
  !! Loop.reorder ~order:((add_prefix "c" dims) @ (add_prefix "b" dims) @ idims) [step; cFor "cX"];

  bigstep "Introduce nbThreads and idThread";
  !! Sequence.insert (expr "#include \"omp.h\"") [tFirst; dRoot];
  !! Variable.insert ~const:false ~name:"nbThreads" ~typ:(atyp "int") ~value:(lit "4") [tBefore; cVarDef "nbCells"];
  (* !! Omp.declare_num_threads "nbThreads"; *) (* TEMPORARY *)
  (* !! Omp.get_num_threads "nbThreads" [tFirst; cTopFunDef "main"; dBody]; *)
  !! Omp.get_thread_num "idThread" [tBefore; step; cFor "iX"];
  !! Trace.reparse();

  bigstep "Parallelize the code using concurrent operations (they are subsequently eliminated)";
  (* BEAUTIFY: should be integrated earlier with "Decompose the loop to allow for parallelization per blocks";*)
  !! Omp.atomic None [tBefore; step; cLabel "charge"; cWrite ()]; (* BEAUTIFY: Instr.set_atomic, and use cOR *)
  !! Expr.replace_fun "bag_push_concurrent" [step; cFun "bag_push"];
  !! Omp.parallel_for ~clause:[Collapse 3] [tBefore; cFor "bX"];
  !! Omp.parallel_for [tBefore; step; cFor "idCell" ~body:[cVar "sum"]];

  (* Part 3: refinement of the parallelization to eliminate atomic operations *)

  bigstep "Duplicate the charge of a corner for each of the threads";
  !! Matrix.delocalize "depositCorners" ~into:"depositThreadCorners" ~indices:["idCell"; "idCorner"]
      ~init_zero:true ~dim:(expr "nbThreads") ~index:"idThread" ~acc_in_place:true ~ops:delocalize_sum ~use:(Some (expr "idThread"))
      [cLabel "core"];
      (* BEAUTIFY:
      ~decl_target:[tAfter; dRoot; cStrict; cVarDef "deposit"]
      ~alloc_target:(Some [tAfter; cTopFunDef "allocateStructures"; cWriteVar "deposit"])
      ~dealloc_target:(Some[tAfter; cTopFunDef "allocateStructures"; cWriteVar "deposit"])*)
  !! Instr.delete [cFor "idCell" ~body:[cCellWrite ~base:[cVar "depositCorners"] ~rhs:[cDouble 0.] ()]];
  !! Variable.init_detach [nbMulti; step; cVarDef ~regexp:true "deposit.*Corners"];
  !! Instr.move_out ~dest:[tAfter; cVarDef "deposit"] [nbMulti; step; cVarDef ~regexp:true "deposit.*Corners"];
  !! Instr.move_out ~dest:[tAfter; cTopFunDef "allocateStructures"; cWriteVar "deposit"] [nbMulti; cWriteVar ~regexp:true "deposit.*Corners"];
  !! Instr.move_out ~dest:[tAfter; cTopFunDef "deallocateStructures"; cFun "free" ~args:[[cVar "field"]]] [nbMulti; step; cFun "MFREE"];
    (* BEAUTIFY:  ~dest:[tLast; cTopFunDef "allocateStructures"; dBody]  *)
  (* TODO: mark the accumulation loop as "depositSum" at the delocalize step
      !! Omp.parallel_for [tBefore; step; cMark "depositSum"]]; *)
  !! Instr.delete [step; cLabel "charge"; cOmp()]; (* BEAUTIFY: Instr.set_nonatomic ; also cPragma is needed *)

  bigstep "Introduce private and shared bags, and use shared ones only for particles moving more than one cell away";
  !! Matrix.delocalize "bagsNext" ~into:"bagsNexts" ~dim:(lit "2") ~indices:["idCell"] ~last:true
    ~alloc_instr:[cFunDef "allocateStructures"; cWriteVar "bagsNext"]
    ~index:"bagsKind" ~ops:delocalize_bag [cLabel "core"];
  !! Variable.insert_list_same_type (atyp "const int") [("PRIVATE", lit "0"); ("SHARED", lit "1")] [tFirst; step; dBody];
  !! Instr.delete [step; cFor "idCell" ~body:[cFun "bag_swap"]];
  !! Variable.exchange "bagsNext" "bagsCur" [nbMulti; step; cFor "idCell"];
  !! Instr.move_out ~dest:[tBefore; cTopFunDef "step"] [step; cOr [[cVarDef "PRIVATE"]; [cVarDef "SHARED"]]];
  !! Instr.delete [cOr[[cVarDef "bagsNext"];[cWriteVar "bagsNext"];[cFun ~regexp:true "\\(free\\|bag.*\\)" ~args:[[cVar "bagsNext"]]]]];
  !! Variable.insert ~typ:(atyp "coord") ~name:"co" ~value:(expr "coordOfCell(idCell2)") [tAfter; step; cVarDef "idCell2"];
  let pushop = cFun "bag_push_concurrent" in
  let force_concurrent_push = false in (* TEMPORARY *)
  let push_cond = if force_concurrent_push then trm_lit (Lit_bool false) else trm_ands (map_dims (fun d ->
         expr ~vars:[d] "(co.i${0} >= b${0} - halfBlock && co.i${0} < b${0} + block + halfBlock)
                      || (b${0} == 0 && co.i${0} >= grid${0} - halfBlock)
                      || (b${0} == grid${0} - block && co.i${0} < halfBlock)")) in
  !! Variable.insert ~typ:(atyp "bool") ~name:"isDistFromBlockLessThanHalfABlock"
      ~value:push_cond [tBefore; step; pushop];
  !! Flow.insert_if ~cond:(var "isDistFromBlockLessThanHalfABlock") [step; pushop];
  !! Specialize.any (expr "PRIVATE") [step; cIf(); dThen; pushop; cAny];
  !! Specialize.any (expr "SHARED") [step; cIf(); dElse; pushop; cAny];
  !! Expr.replace_fun "bag_push_serial" [step; cIf(); dThen; pushop];
  !! Trace.reparse ();
  !! Variable.init_detach [nbMulti; step; cVarDef "bagsNexts"]; (*BEAUTIFY: attach to the delocalize op *)
  !! Instr.move_out ~dest:[tAfter; cVarDef "deposit"] [nbMulti; step; cVarDef "bagsNexts"];
  !! Instr.move_out ~dest:[tAfter; cTopFunDef "allocateStructures"; cWriteVar "deposit"] [nbMulti; cWriteVar "bagsNexts"];
  !! Instr.move_out ~dest:[tAfter; cTopFunDef "deallocateStructures"; cFun "free" ~args:[[cVar "field"]]] [nbMulti; step; cFun "MFREE"];
  !! Instr.move_out ~dest:[tAfter; cTopFunDef "allocateStructures"; cFor ""] [nbMulti;step; cFor "idCell" ~body:[cFun "bag_init_initial"]];
  !! Instr.move_out ~dest:[tBefore; cTopFunDef "deallocateStructures"; cFor ""] [nbMulti;step; cFor "idCell" ~body:[cFun "bag_free_initial"]];
  !! Instr.delete [step; cVarDef "contribs"];

  bigstep "Parallelize and optimize loops that process bags";
  !! Loop.fusion ~nb:2 [step; cFor "idCell" ~body:[cFun "bag_append"]];
  !! Omp.parallel_for [tBefore; occIndex 1; step; cFor "idCell"]; (* BEAUTIFY: use label to refer to the loop *)
  !! Function.use_infix_ops ~indepth:true [step; dBody]; (* LATER: move to the end of an earlier bigstep *)

  (* Part 4: Vectorization *)

  bigstep "Loop splitting: process speeds, process positions, deposit particle and its charge";
  (* Unrolling coeff computation loops and inlining coeff writes *)
  !! Loop.unroll [occFirst; step; cFor "i"; cFor "k"];
  !! Loop.unroll [stepLF; cFor "i"; cFor "k"];
  !! Instr.inline_last_write [nbMulti; steps; cCellRead ~base:[cFieldRead ~base:[cVar "coeffs"] ()] ()];
  !! iter_dims (fun d ->
       Instr.inline_last_write [steps; cCellWrite ~base:[cFieldRead ~field:("itemsSpeed"^d) ()] (); cReadVar ("fieldAtPos"^d)];);
  !! Variable.inline [nbMulti; steps; cVarDef ~regexp:true "fieldAt.*"];
  (* BEAUTIFY: arith simpl should remove the  0. + .. *)
  (* BEAUTIFY: in the printing functions we have things like
        ToBuffer.pretty 0.9 80 b (ast_to_doc t);
     we should replace 80 with (!Flags.code_print_width),
     so that we can customize it, e.g. to 120 for pic_demo.ml *)
  !! Trace.reparse();
  !! Variable.to_nonconst [step; cVarDef "idCell2"];
  !! Loop.hoist ~array_size:(Some (expr "CHUNK_SIZE")) [step; cVarDef "idCell2"];
     let dest = [tBefore; step; cVarDef "isDistFromBlockLessThanHalfABlock"] in
  !! Instr.copy ~dest [step; cVarDef "idCell2"];
  !! Instr.move ~dest [step; cVarDef "co"];
  !! Loop.fission [nbMulti; tBefore; step; cOr [[cVarDef "pX"]; [cVarDef "rX1"]]];
  !! Variable.ref_to_pointer [occFirst; step; cVarDef "idCell2"];
  !! Variable.ref_to_var [occLast; step; cVarDef "idCell2"];

  (* BEAUTIFY
      int *idCell2 = &idCell2_step[i];
      *idCell2 = MINDEX3(gridX, gridY, gridZ, iX2, iY2, iZ2);
    YES THIS one you can replace with
       idCell2_step[i] = MINDEX3(gridX, gridY, gridZ, iX2, iY2, iZ2);
      just by doing the variable inlining (cancellation of *& should be done on the fly
      by variable inline)
  *)

  (* BEAUTIFY:
    (double [star])MMALLOC_ALIGNED1(nbCells, sizeof(double), 64);
    SHOULD BE
    (double[star]) MMALLOC_ALIGNED1(nbCells, sizeof(double), 64);
    there are two spaces to modify *)

  bigstep "Vectorization";
  (* TODO: Clean up  *)
  !! Sequence.insert (expr "#include \"stdalign.h\"") [tFirst; dRoot]; (* BEAUTIFY: Align.header [] *)
  !! Loop.fission [tBefore; occLast; step; cFor "idCell"; cFor "idCorner"; cFor "idThread"];
  !! Loop.swap [occLast; cFor "idCell"; cFor "idCorner" ~body:[cFor "idThread"]];
  !! Omp.simd [nbMulti; tBefore;cFor "idCell" ~body:[cFor "bagsKind"]; cFor "idCorner"];
  !! Omp.simd ~clause:[Aligned (["coefX"; "coefY"; "coefZ"; "signX"; "signY"; "signZ"], align)] [tBefore; step; cLabel "charge"];
  !! Omp.simd [tBefore; stepLF; cFor "i"];
  !! Label.remove [step; cLabel "charge"];
  !! Align.def (lit "64") [nbMulti; cOr [[cStrict; cVarDef ~regexp:true "\\(coef\\|sign\\)."];
                                         [step; cVarDef "idCell2_step"];
                                         [cStrict; cVarDef ~substr:true "deposit"]]];
  !! Struct.align_field (lit "64") ("items.") [cTypDef "chunk"];
  !! List.iter (fun occ -> Omp.simd [occIndex occ; tBefore; step; cFor "i"]) [0;1]; (* BEAUTIFY: occIndices *)
  !! Function.inline [step; cFun "cellOfCoord"];
  !! Align.alloc (lit "64") [nbMulti; cTopFunDef "allocateStructures"; cMalloc ()];
  !! Function.inline [nbMulti; step; cFor "idCorner" ~body:[cVar "depositCorners"]; cMindex ~d:2 ()];
  !! Instr.inline_last_write [step; cFor "idThread"; cReadVar "__TEMP_Optitrust"];
  !! Instr.delete [step; cVarDef "__TEMP_Optitrust"];
  !! Function.inline [nbMulti; step; cFor "idCorner" ~body:[cVar "depositCorners"]; cMindex ~d:3 ()];
  !! Variable.init_attach [step; cVarDef "__TEMP_Optitrust"];
  !! Variable.inline [step; cVarDef "__TEMP_Optitrust"];
  !! Function.inline [nbMulti; step; cFor "k"; cMindex ~d:3 ()];
)
(*


  LATER:
  const int idCell = (iX * gridY + iY) * gridZ + iZ;
  could be
  cellOfCoord(iX, iY, iZ)
  if the user provides the name for this function



??
use unfold in
!! Variable.inline [nbMulti; step; cVarDef ~regexp:true "factor."];
  with a more precise target


  LATER: keep in mind
  Align.header => adds  #include <stdalign.h> to the top of the ast
  Align.assume "t" =>
     insert instruction
     t = __builtin_assume_aligned(t, 64);
  Align.def 64 [cVarDef "x"]
    => add the attribute to the type of the definition

  LATER
  Flags.print_coumpound_expressions



 TODO:
  introduce
    int id = *idCell2;
  just before coordOfCell( *idCell2)
    should fold3 occurences.

    inline idCell2  everywhere


*)
(* LATER !! Function.beta ~indepth:true [dRoot]; try in a unit test with two beta reductions to do *)

(* TODO: res should be field_at_corners *)



(* if not doublepos then begin
    bigstep "Turn positions into floats";
    !! Cast.insert (atyp "float") [sExprRegexp ~substr:true "p.2 - i.2"];
    LATER: target the [sRexegp "c->itemsPos[[.]] = ")  or iter_dims  or   ==>best: cOr (map_dims .. )
    !! Struct.update_fields_type "itemsPos." (atyp "float") [cTypDef "chunk"];
    LATER: type particle would need to be converted too
       const vect pos = {x, y, z};
       would need cast around the values
  end;
*)


  (* LATER: bigstep "Introduce matrix operations, and prepare loop on charge deposit";
   might be useful to group this next to the reveal of x/y/z *)

  (* LATER: menhir should support "res[]" syntax *)

  (* LATER: !! Expr.replace ~reparse:false (expr "MINDEX2(nbCells, 8, idCell2, k)")
           [step; cLabel "charge"; cFun "mybij"];
      instead use: sExpr "mybij(nbCorners, nbCells, indicesOfCorners(idCell2).v[k], k)"
       ARTHUR: simplify mybij calls in the sum *)


  (* LATER !! Instr.delete [cOr[[cVarDef "bagsNext"];[ cKindInstr; cVar "bagsNext"]]]; *)
  (* LATER !! Variable.delete [cVarDef "bagsNext"] ==> shorthand for above *)


  (* cleanup *)
  (* TODO: delocalize_obj ~labels:["allocBagsNext","","deallocBagsNext"] => creates a sequence with malloc and init loop *)
  (* TODO: move allocBagsNext and deallocBagsNext labelled blocks *)



  (* LATER: fission should automatically do the duplication of references when necessary *)

  (*  BEAUTIFY: the Omp operations should  add a tBefore themselves
    !! Omp.simd [] [tBefore; step;cFor "i"]; *)(* TODO: Fix the issue with the last loop *)


(* TODO: Fix the bug with insertion of variables when using tBefore and tFirst *)


(* !! Variable.insert_list ~reparse:true ~defs:(List.rev ( TODO
         ["const double", "factorC", expr "particleCharge * stepDuration * stepDuration / particleMass"]
       @ (map_dims (fun d -> "const double", ("factor" ^ d), expr ("factorC / cell" ^ d)))))
     [tFirst; step; dBody]; *)

(* LATER: rename pic_demo.c to pic_naive.c *)

(*LATER halfBlock=expr "block/2"*)

(* TODO

  struct vect {
      alignas(64) double x;   pattern="x"     newtypetouse=(typ_update "double")
      alignas(64) int y;
      }int : typ) : Target.Transfo.t =
  applyto_fields_type pattern (fun _ -> ty)

let Ast.typ_alignas (align:int) (ty : typ) =
  same typ with attributes typ_alignas added

let Combi_Struct.align_field (align:int) (pattern : string) =
  Struct.applytofields_type (fun ty -> typ_alignas align ty)

*)


(* LATER

// after createParticle, add applyScalingShifting(true)
// after cFor "idStep" in main, add applyScalingShifting(false)
void applyScalingShifting(bool dir) { // dir=true at entry, dir=false at exit
 for (int idCell = 0; idCell < nbCells; idCell++) {
    bag* b = &bagsCur[idCell];
    bag_iter bag_it;
    for (particle* p = bag_iter_begin(&bag_it, b); p != NULL; p = bag_iter_next_common(&bag_it, false)) {
      p->pos.x = p->pos.x;
      p->pos.y = p->pos.y;
      p->pos.z = p->pos.z;
      p->speed.x = p->speed.x;
      p->speed.y = p->speed.y;
      p->speed.z = p->speed.z;
    }
 }
}
/*
      p->pos.x = (p->pos.x + ix) * cellX;
      p->pos.y = p->pos.y;
      p->pos.z = p->pos.z;
      p->speed.x = p->speed.x * cellX / stepDuration;
      p->speed.y = p->speed.y;
      p->speed.z = p->speed.z;
*/

*)

(* LATER
  !! List.iter (fun occ -> Omp.parallel_for [occIndex occ; tBefore; step; cFor "idCell"]) [0(*;2*)]; (* BEAUTIFY: occIndices *)
  *)



(* BEAUTIFY implement and use: insert_list_same_type *)

(* BEAUTIFY
  missing spaces in:
    return (vect){d * v.x, d * v.y, d * v.z};
  should be
   return (vect) { d * v.x, d * v.y, d * v.z };

  too many spaces in:
    const  double cellVolume = cellX * cellY * cellZ;
  probably due to the way you print aliasas attributes
   *)


(* LATER: keep this code, it might be useful in the future
     let wrapPow_def = "int wrapPowersOfTwo(int gridSize, int a) {return a & (gridSize - 1);}" in
  !! Sequence.insert ~reparse:true (stmt wrapPow_def) [tBefore; step];
  !! Function.inline [nbMulti; step; cFun "wrapPowersOfTwo"];*)


(* LATER: the C standard parses
      x & y - 1
  as
     x & (y - 1)
  but this is very confusing, so we should always put parentheses around nontrivial arguments of & and | operators.
  *)


