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
let delocalize_bag = Local_obj ("bag_init", "bag_append", "bag_free")
let align = 64

(* Grab the "usechecker" flag from the command line *)
let usechecker = ref false
let _= Run.process_cmdline_args
  [("-usechecker", Arg.Set usechecker, " use -DCHECKER as preprocessor flag")]
  (* LATER: use a generic -D flag for optitrust *)
let usechecker = !usechecker

(* let _ = Printf.printf "usechecker=%d\n" (if usechecker then 1 else 0) *)

(* UNCOMMENT THE LINE BELOW FOR WORKING ON THE VERSION WITH THE CHECKER *)
(* let usechecker = true  *)

let grid_dims_power_of_2 = true

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

(* let _ = Flags.code_print_width := 120 *)

let _ = Run.script_cpp ~parser:Parsers.Menhir ~prepro ~inline:["pic_demo.h";"bag.hc";"particle.hc";"optitrust.h";"bag_atomics.h";"bag.h-"] (fun () ->

  Printf.printf "CHECKER=%d\n" (if usechecker then 1 else 0);

  (* Part 1: sequential optimizations *)

  bigstep "Optimization and inlining of [matrix_vect_mul]";
  let ctx = cTopFunDef "matrix_vect_mul" in
  !! Function.inline [ctx; cFuns ["vect_mul"; "vect_add"]];
  !! Struct.set_explicit [nbMulti; ctx; cWriteVar "res"];
  !! Loop.fission ~split_between:true [ctx; cFor "idCorner"];
  !! Loop.unroll [nbMulti; ctx; cFor "idCorner"];
  !! Instr.accumulate ~nb:8 [nbMulti; ctx; sInstrRegexp ~substr:true "res.*\\[0\\]"];
  !! Function.inline ~delete:true [nbMulti;cFun "matrix_vect_mul"];

  bigstep "Optimization in [cornerInterpolationCoeff], before it is inlined";
  let ctx = cTopFunDef "cornerInterpolationCoeff" in
  let ctx_rv = cChain [ctx; sInstr "r.v"] in
  !! Rewrite.equiv_at "double a; ==> 1. - a == (1. + -1. * a)" [nbMulti; ctx; cVarDef ~regexp:true "c."; cInit()];
  !! Rewrite.equiv_at "double a; ==> a == (0. + 1. * a)" [nbMulti; ctx_rv; cVar ~regexp:true "r."];
  !! Variable.inline [nbMulti; ctx; cVarDef ~regexp:true "c."];
  !! Variable.intro_pattern_array ~const:true ~pattern_aux_vars:"double rX, rY, rZ"
      ~pattern_vars:"double coefX, signX, coefY, signY, coefZ, signZ"
      ~pattern:"(coefX + signX * rX) * (coefY + signY * rY) * (coefZ + signZ * rZ)"
      [nbMulti; ctx_rv; dRHS];
  !! Instr.move_out ~dest:[tBefore; ctx] [nbMulti; ctx; cVarDef ~regexp:true "\\(coef\\|sign\\)."];
  !! Loop.fold_instrs ~index:"idCorner" [cTopFunDef "cornerInterpolationCoeff"; cCellWrite ~base:[cVar "r"] ()];

  bigstep "Eliminate an intermediate storage by reusing an existing one";
  !! Variable.reuse ~space:(expr "p->speed") [step; cVarDef "speed2" ];
  !! Variable.reuse ~space:(expr "p->pos") [step; cVarDef "pos2"];

  bigstep "Reveal write operations involved in the manipulation of particles and vectors";
  let ctx = cTopFunDefs ["bag_push_serial";"bag_push_concurrent"] in
  !! List.iter (fun typ -> Struct.set_explicit [nbMulti; ctx; cWrite ~typ ()]) ["particle"; "vect"];
  !! Function.inline [steps; cFuns ["vect_mul"; "vect_add"]];
  !! Trace.reparse ();
  !! Struct.set_explicit [nbMulti; stepsReal; cFieldWrite ~base:[cVar "p"] ()];
  !! Function.inline ~delete:true ~vars:(AddSuffix "${occ}") [nbMulti; step; cFun "wrapArea"];
  !! Variable.inline [nbMulti; step; cVarDef ~regexp:true "[xyz]."];

  bigstep "Inlining of [cornerInterpolationCoeff] and [accumulateChargeAtCorners]";
  !! Function.inline [nbMulti; cTopFunDef "cornerInterpolationCoeff"; cFun ~regexp:true "relativePos."];
  !! Function.inline [step; cFun "accumulateChargeAtCorners"];
  !! Function.inline ~vars:(AddSuffix "2") [step; cFun "idCellOfPos"];
  !! List.iter (fun f -> Function.inline ~vars:(AddSuffix "${occ}") [nbMulti; f; cFun "cornerInterpolationCoeff"])
     stepsl;
  !! iter_dims (fun d -> Variable.reuse ~space:(var ("i" ^ d ^ "2")) [step; cVarDef ("i" ^ d ^ "1")]);
  !! Trace.reparse();

  bigstep "Simplification of the deposit of charge";
  !! Sequence.intro ~mark:"fuse" ~start:[step; cVarDef "contribs"] ();
  !! Loop.fusion_targets [cMark "fuse"];
  !! Instr.inline_last_write [step; cCellRead ~base:[cFieldRead ~base:[cVar "contribs"] ()] ()];

  bigstep "Low level iteration on chunks of particles";
  !! Function.inline [steps; cFuns ["bag_iter_begin"; "bag_iter_destructive_begin"]];
  !! Sequence.intro_on_instr [steps; cFor_c ""; dBody];
  !! Function.uninline ~fct:[cTopFunDef "bag_iter_ho_basic"] [steps; cVarDef "bag_it"];
  !! Expr.replace_fun ~inline:true "bag_iter_ho_chunk" [steps; cFun "bag_iter_ho_basic"];
  !! List.iter (fun f -> Function.beta ~indepth:true [f]) stepFuns;
  !! Instr.delete [nbMulti; cTopFunDefAndDecl ~regexp:true "bag_iter.*"];

  bigstep "Elimination of the pointer on a particle, to prepare for aos-to-soa";
  !! Instr.inline_last_write [nbMulti; steps; cRead ~addr:[cStrictNew; cVar "p"] ()];

  bigstep "Preparation for AOS-TO-SOA";
  !! Struct.set_explicit [step; cVarDef "p2"];
  !! Struct.set_explicit [nbMulti; step; cFieldWrite ~base:[cVar "p2"] ~regexp:true ~field:"\\(speed\\|pos\\)" ()];

  bigstep "AOS-TO-SOA";
  !! List.iter (fun f -> Struct.reveal f [cTypDef "particle"]) ["speed"; "pos"];
  !! Struct.reveal "items" [cTypDef "chunk"];

  bigstep "Apply scaling factors on the electric field";
  !! Struct.to_variables [steps; cVarDef "fieldAtPos"];
  !! Variable.insert_list_same_type ~reparse:true (ty "const double") (["factorC", expr "particleCharge * stepDuration * stepDuration / particleMass"]
      @ (map_dims (fun d -> ("factor" ^ d, expr ("factorC / cell" ^ d))))) [tBefore; steps; cFor "idCell" ~body:[cFor "i"]];
  !! Function.inline ~delete:true [steps; cFun "getFieldAtCorners"];
  !! Variable.rename ~into:"field_at_corners" [step; cVarDef "res"];
  !! Struct.set_explicit [steps; cFor "idCorner"; cCellWrite ~base:[cFieldRead ~base:[cVar "field_at_corners"] ()] ()];
  !! iter_dims (fun d ->
      let d1 = String.lowercase_ascii d in
      Accesses.scale ~factor:(var ("factor" ^ d)) [steps; cFor "idCorner"; cFieldWrite ~field:d1 ()]);
  !! iter_dims (fun d ->
       Accesses.scale ~factor:(var ("factor" ^ d)) [steps; cVarDef "accel"; cReadVar ("fieldAtPos" ^ d)]);
  !! Variable.unfold [step; cVarDef  "factorC"];
  !! Variable.unfold ~at:[cVarDef "accel"] [nbMulti; step; cVarDef ~regexp:true "factor."];
  !! Arith.(simpl ~indepth:true expand) [nbMulti; steps; cVarDef "accel"];

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
      Accesses.scale ~factor:(var ~mut:true ("cell"^d)) [addPart; cFieldRead ~field:d_lc ~base:[cVar "pos"] ()]);
  !! iter_dims (fun d ->
     Accesses.scale ~neg:true ~factor:(var ~mut:true ("cell"^d))
         [nbMulti; steps; cOr [[sExprRegexp ("c->itemsPos" ^ d ^ "\\[i\\]")]; [cFieldWrite ~field:("pos"^d)()]]]);

  bigstep "Simplify arithmetic expressions after scaling";
  !! Trace.reparse();
  !! Variable.inline [steps; cVarDef "accel"];
  !! Arith.with_nosimpl [nbMulti; steps; cFor "idCorner"] (fun () ->
       Arith.(simpl ~indepth:true expand) [nbMulti; steps]);
  (* !! Function.use_infix_ops ~indepth:true [step; dBody]; *)

      (* BEAUTIFY: remind me why we can't do a infixop just here? it would be very pretty;
          (even if we have to undo it later) *)

  bigstep "Enumerate grid cells by coordinates";
  !! Instr.read_last_write ~write:[cWriteVar "nbCells"] [step; cFor "idCell" ~body:[cFor "i"]; cReadVar "nbCells"];
  !! Loop.grid_enumerate ~indices:(map_dims (fun d -> "i"^d)) [step; cFor "idCell" ~body:[cFor "idCorner"]];

  bigstep "Code cleanup in preparation for shifting of positions";
  !! iter_dims (fun d ->
      Variable.bind ~const:true ~typ:(Some (ty "double")) ("p" ^ d ^ "2") [occLast;step; cCellWrite ~base:[cFieldRead ~field:("itemsPos" ^ d) ()] (); dRHS];
      Variable.bind ~const:true ("p" ^ d ) ~typ:(Some (ty "double")) [occFirst;step; cCellWrite ~base:[cFieldRead ~field:("itemsPos" ^ d) ()] (); dRHS]);
  !! iter_dims (fun d ->
      Instr.inline_last_write [step; cVarDef ~regexp:true "p.2"; cCellRead ~base:[cFieldRead ~field:("itemsPos" ^ d) ()]()]);
  !! iter_dims (fun d ->
      Instr.read_last_write [step; cVarDef ~regexp:true "i.2"; cCellRead ~base:[cFieldRead ~field:("itemsPos" ^ d) ()]()]);
  !! Instr.(gather_targets ~dest:GatherAtFirst) [step; cVarDef ~regexp:true ("\\(i.*2\\|p.2\\|i.2\\)")];
  !! Instr.(gather_targets ~dest:(GatherAt [tBefore; cVarDef "p2"])) [step; cVarDef ~regexp:true "r.1"];

  bigstep "Shifting of positions: make positions relative to the containing cell";
  !! Instr.move ~dest:[tBefore; addPart; cVarDef "p"] [addPart; cVarDef "idCell"];
  !! List.iter (fun tg ->
      Variable.insert ~typ:(ty "coord") ~name:"co" ~value:(expr "coordOfCell(idCell)") tg)
      ([ [tAfter; addPart; cVarDef "idCell"] ] @ onlychecker [tFirst; repPart; cFor "idCell"; dBody]); (* LATER: make cOr work for targetBetweens (hard) *)
  !! iter_dims (fun d ->
      Accesses.shift ~neg:true ~factor:(expr ("co.i"^d)) [cOr (
        [[addPart; cFieldWrite ~field:("pos"^d) ()]] @
        (onlychecker [repPart; cCellRead ~base:[cFieldRead ~field:("itemsPos" ^ d) ()] ()]) )] );

  !! iter_dims (fun d ->
      Accesses.shift ~neg:true ~factor:(var ("i" ^ d ^ "0")) [stepsReal; cVarDef ~regexp:true "r.0"; cCellRead ~base:[cFieldRead ~field:("itemsPos" ^ d) ()]()];
      Accesses.shift ~neg:true ~factor:(var ("i" ^ d)) [step; cVarDef ~regexp:true ("p" ^ d); cCellRead ~base:[cFieldRead ~field:("itemsPos" ^ d) ()]()];
      Accesses.shift ~neg:true ~factor:(var ("i" ^ d ^ "2")) [step; cCellWrite ~base:[cFieldRead ~field:("itemsPos" ^ d) ()]()];
      Accesses.shift ~neg:true ~factor:(var ("i" ^ d ^ "2")) [step; cVarDef ("r" ^ d ^ "1"); cCellRead ~base:[cFieldRead ~field:("itemsPos" ^ d) ()]()]);

  bigstep "Simplify arithmetic expressions after shifting of positions";
  !! Rewrite.equiv_at ~glob_defs:"double fwrap(double, double);\n" "double x, y, z; ==> (fwrap(x,y)/z) == (fwrap(x/z, y/z))" [cVarDef ~regexp:true "p.2"; cInit()];
  !! iter_dims (fun d ->
      Instr.read_last_write ~write:[cTopFunDef "computeConstants"; cWriteVar ("cell"^d)] [nbMulti;step; cFun "fwrap";cReadVar ("cell"^d)];);

  !! Arith.with_nosimpl [nbMulti; stepsReal; cFor "idCorner"] (fun () ->
       Arith.(simpl ~indepth:true expand) [nbMulti; stepsReal]);
  !! Variable.inline [stepsReal; cVarDef ~regexp:true "r.."];
  !! Instr.delete [nbMulti; stepsReal; cVarDef ~regexp:true "i.0"];

  if doublepos then begin
    bigstep "Turn positions into floats, decreasing precision but allowing to fit a larger number of particles";
    !! Cast.insert (ty "float") [sExprRegexp ~substr:true "p.2 - i.2"];
    !! Struct.update_fields_type "itemsPos." (ty "float") [cTypDef "chunk"];
  end;

  bigstep "Replacement of the floating-point wrap-around operation with an integer wrap-around";
   let fwrapInt = "double fwrapInt(int m, double v) {
      const int q = int_of_double(v);
      const double r = v - q;
      const int j = wrap(m, q);
      return j + r;
    }" in
  !! Function.insert ~reparse:true fwrapInt [tBefore; step];
  !! Expr.replace_fun "fwrapInt" [nbMulti; step; cFun "fwrap"];
  !! iter_dims (fun d ->
      Function.inline ~vars:(AddSuffix d) [step; cVarDef ("p"^d^"2"); cFun "fwrapInt"]);
  if grid_dims_power_of_2 then
    !! Rewrite.equiv_at ~glob_defs:"int wrap(int, int);\n" "int a, b; ==> wrap(a,b) == (b & (a -1))" [nbMulti; step; cFun "wrap"];

  bigstep "Simplification of computations for positions and destination cell";
  !! iter_dims (fun d -> Expr_basic.replace (var ("j"^d)) [step; cVarDef ("i"^d^"2");cInit()];);
  !! Variable.inline_and_rename [nbMulti; step; cVarDef ~regexp:true "i.2" ];
  !! Variable.inline [nbMulti; step; cVarDef ~regexp:true "p.2"];
  !! Arith.(simpl_rec expand) [nbMulti; step; cCellWrite ~base:[cFieldRead ~regexp:true ~field:("itemsPos.") ()] ()];

  bigstep "Introduce matrix operations, and prepare loop on charge deposit";
  !! Label.add "core" [step; cFor "iX" ];
  !! Matrix.intro_mops (var ~mut:true "nbCells") [nbMulti; cVarDef ~regexp:true "\\(deposit\\|bagsNext\\)"];
  !! Label.add "charge" [step; cFor "idCorner" ~body:[cVar "deposit"]];
  !! Variable.inline [occLast; step; cVarDef "indices"];

  bigstep "Duplicate the charge of a corner for the 8 surrounding cells";
  let alloc_instr = [cTopFunDef "allocateStructures"; cWriteVar "deposit"] in
  !! Matrix.delocalize "deposit" ~into:"depositCorners" ~last:true ~indices:["idCell"] ~init_zero:true
     ~labels:["alloc"; ""; "dealloc"] ~dealloc_tg:(Some [cTopFunDef ~regexp:true "dealloc.*"; cFor ""])
     ~dim:(var "nbCorners") ~index:"idCorner" ~acc:"sum" ~ops:delocalize_sum ~use:(Some (var "idCorner")) ~alloc_instr [cLabel "core"];

  bigstep "Apply a bijection on the array storing charge to vectorize charge deposit";
  let bij_def =
      "int bij(int nbCells, int nbCorners, int idCell, int idCorner) {
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
  !! Function.insert bij_def [tBefore; step];
  !! Matrix.biject "bij" [cVarDef "depositCorners"];
  !! Expr.replace ~reparse:false (expr "MINDEX2(nbCells, nbCorners, idCell2, idCorner)")
      [step; sExpr "bij(nbCells, nbCorners, indicesOfCorners(idCell2).v[idCorner], idCorner)"];

  (* Part 2: parallelization *)

  bigstep "Decompose the loop to allow for parallelization per blocks";
  !! Variable.insert_list_same_type (ty "const int") [("block", lit "2"); ("halfBlock", (lit "1"))] [tBefore; cVarDef "nbCells"];
  !! iter_dims (fun d -> let bd = "b"^d in
      Loop.tile (var ~mut:true "block") ~bound:TileBoundDivides ~index:("b"^d) [step; cFor ("i"^d)];
      Loop.color (lit "2") ~index:("c"^d) [step; cFor bd] );
  !! Loop.reorder ~order:((add_prefix "c" dims) @ (add_prefix "b" dims) @ idims) [step; cFor "cX"];
  !! Expr.replace_fun "bag_push_concurrent" [step; cFun "bag_push"];
  !! Instr.set_atomic [step; cLabel "charge"; cWrite ()];
  !! Omp.parallel_for ~collapse:3 [step; cFor "bX"];

  bigstep "Introduce nbThreads and idThread";
  !! Omp.header ();
  !! Variable.insert ~const:false ~name:"nbThreads" ~typ:(ty "int") [tBefore; cVarDef "nbCells"];
  !! Omp.get_thread_num "idThread" [step; cFor "iX"];
  !! Omp.set_num_threads ("nbThreads") [tFirst; cTopFunDef "main"; dBody];
  !! Trace.reparse();

  bigstep "Duplicate the charge of a corner for each of the threads";
  let alloc_instr = [cTopFunDef "allocateStructures"; cWriteVar "depositCorners"] in
  !! Matrix.delocalize "depositCorners" ~into:"depositThreadCorners" ~indices:["idCell"; "idCorner"]
      ~init_zero:true ~dim:(var ~mut:true "nbThreads") ~index:"idThread" ~acc_in_place:true ~ops:delocalize_sum ~use:(Some (var "idThread"))
      ~labels:["alloc"; ""; "dealloc"] ~alloc_instr ~dealloc_tg:(Some [cTopFunDef ~regexp:true "dealloc.*"; cFor ""])
      [cLabel "core"];
  !! Instr.delete [cFor "idCell" ~body:[cCellWrite ~base:[cVar "depositCorners"] ~rhs:[cDouble 0.] ()]];
  !! Instr.delete [step; cLabel "charge"; cOmp()]; (* BEAUTIFY: Instr.set_nonatomic ; also cPragma is needed *)

  bigstep "Introduce private and shared bags, and use shared ones only for particles moving more than one cell away";
  !! Matrix.delocalize "bagsNext" ~into:"bagsNexts" ~dim:(lit "2") ~indices:["idCell"] ~last:true
    ~alloc_instr:[cTopFunDef "allocateStructures"; cWriteVar "bagsNext"]
    ~labels:["alloc"; ""; "dealloc"] ~dealloc_tg:(Some [cTopFunDef ~regexp:true "dealloc.*"; cFor ""])
    ~index:"bagsKind" ~ops:delocalize_bag [cLabel "core"];
  !! Variable.insert_list_same_type (ty "const int") [("PRIVATE", lit "0"); ("SHARED", lit "1")] [tBefore; step];
  !! Instr.delete [step; cFor "idCell" ~body:[cFun "bag_swap"]];
  !! Variable.exchange "bagsNext" "bagsCur" [nbMulti; step; cFor "idCell"];
  !! Instr.delete [cOr[[cVarDef "bagsNext"];[cWriteVar "bagsNext"];[cFun ~regexp:true "\\(free\\|bag.*\\)" ~args:[[cVar "bagsNext"]]]]];
  !! Variable.insert ~typ:(ty "coord") ~name:"co" ~value:(expr "coordOfCell(idCell2)") [tAfter; step; cVarDef "idCell2"];
  let pushop = cFun "bag_push_concurrent" in
  let force_concurrent_push = false in (* TEMPORARY *)
  let push_cond = if force_concurrent_push then trm_lit (Lit_bool false) else trm_ands (map_dims (fun d ->
         expr ~vars:[d] "(co.i${0} >= b${0} - halfBlock && co.i${0} < b${0} + block + halfBlock)
                      || (b${0} == 0 && co.i${0} >= grid${0} - halfBlock)
                      || (b${0} == grid${0} - block && co.i${0} < halfBlock)")) in
  !! Variable.insert ~typ:(ty "bool") ~name:"isDistFromBlockLessThanHalfABlock"
      ~value:push_cond [tBefore; step; pushop];
  !! Flow.insert_if ~cond:(var "isDistFromBlockLessThanHalfABlock") [step; pushop];
  !! Specialize.any (var ~mut:true "PRIVATE") [step; cIf(); dThen; pushop; cAny];
  !! Specialize.any (var ~mut:true "SHARED") [step; cIf(); dElse; pushop; cAny];
  !! Expr.replace_fun "bag_push_serial" [step; cIf(); dThen; pushop];

  bigstep "Parallelize and optimize loops that process bags";
  !! Trace.reparse();
  !! Loop.fusion ~nb:2 [step; cFor "idCell" ~body:[cFun "bag_append"]];
  !! Omp.parallel_for [nbMulti;stepsReal; cFor "idCell"];
  !! Function.use_infix_ops ~indepth:true [step; dBody]; (* LATER: move to the end of an earlier bigstep *)

  (* Part 4: Vectorization *)

  bigstep "Loop splitting: process speeds, process positions, deposit particle and its charge";
  (* Unrolling coeff computation loops and inlining coeff writes *)
  !! Loop.unroll [occFirst; step; cFor "i"; cFor "idCorner"];
  !! Loop.unroll [stepLF; cFor "i"; cFor "idCorner"];
  !! Instr.inline_last_write [nbMulti; steps; cCellRead ~base:[cFieldRead ~base:[cVar "coeffs"] ()] ()];
  !! iter_dims (fun d ->
       Instr.inline_last_write [steps; cCellWrite ~base:[cFieldRead ~field:("itemsSpeed"^d) ()] (); cReadVar ("fieldAtPos"^d)];);
  !! Variable.inline [nbMulti; steps; cVarDef ~regexp:true "fieldAt.*"];
  !! Arith.with_nosimpl [nbMulti; stepsReal; cFor "idCorner"] (fun () ->
       Arith.(simpl ~indepth:true expand) [nbMulti; stepsReal; cFor "i"]);
  !! Variable.to_nonconst [step; cVarDef "idCell2"];
  !! Loop.hoist ~array_size:(Some (var "CHUNK_SIZE")) [step; cVarDef "idCell2"];
     let dest = [tBefore; step; cVarDef "isDistFromBlockLessThanHalfABlock"] in
  !! Instr.copy ~dest [step; cVarDef "idCell2"];
  !! Instr.move ~dest [step; cVarDef "co"];
  !! Loop.fission [nbMulti; tBefore; step; cOr [[cVarDef "pX"]; [cVarDef "p2"]]];
  !! Variable.inline [nbMulti; step; cVarDef "idCell2"];
  !! Instr.delete [step; cVarDef "contribs"];


  bigstep "Data alignment";
  !! Align.def (lit "64") [nbMulti; cOr [[cStrict; cVarDef ~regexp:true "\\(coef\\|sign\\)."];
                                         [step; cVarDef "idCell2_step"];
                                         [cStrict; cVarDef ~substr:true "deposit"]]];
  !! Struct.align_field (lit "64") ("items.") [cTypDef "chunk"];
  !! Function.inline [step; cFun "cellOfCoord"];
  !! Align.alloc (lit "64") [nbMulti; cTopFunDef "allocateStructures"; cMalloc ()];

  bigstep "Function inlining on loops to be vectorized";
  let ctx = cFor "idCorner" ~body:[cVar "depositCorners"] in
  !! Function.bind_intro ~fresh_name:"temp_var_${occ}" [nbMulti; ctx; cMindex ()];
  !! Function.inline [nbMulti; ctx; cMindex ()];
  !! Variable.inline [nbMulti; ctx; cVarDef ~regexp:true "temp_var_."];

  bigstep "Vectorization";
  !! Align.header ();
  !! Loop.fission [tBefore; occLast; step; cFor "idCell"; cFor "idCorner"; cFor "idThread"];
  !! Loop.swap [occLast; cFor "idCell"; cFor "idCorner" ~body:[cFor "idThread"]];
  !! Omp.simd [nbMulti; cOr [
    [cFor "idCell" ~body:[cFor "bagsKind"]; cFor "idCorner"];
    [stepLF; cFor "i"];
    [cDiff [[step; cFor "i"]] [[step; cFor "i"  ~body:[cFor "idCorner"]]] ]; (* A / B logic *)
    ]];
  !! Omp.simd ~clause:[Aligned (["coefX"; "coefY"; "coefZ"; "signX"; "signY"; "signZ"], align)] [step; cLabel "charge"];
  !! Label.remove [step; cLabel "charge"];

)
