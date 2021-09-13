open Optitrust
open Target

let _ = Run.script_cpp (fun () ->

  (* One (ExpectedOne is the default) *)
  show [ cFor "j" ];

  (* Multi *)
  show [ nbMulti; cFor "i" ];
  show [ nbMulti; cCall "f" ];
  show [ nbMulti; cFunDef "main"; cFor "i" ];

  (* Nb *)
  show [ nbExact 0; cFunDef "main"; cFor "j" ]; (* zero match *)
  show [ nbExact 1; cFunDef "main"; cFor "i" ];
  show [ nbExact 2; cFor "i" ];

  (* Any *)
  show [ nbAny; cFunDef "main"; cFor "j" ]; (* zero match *)
  show [ nbAny; cFunDef "main"; cFor "i" ];
  show [ nbAny; cFor "i" ];
  show [nbMulti;cAnd [[cReturn];[cTypDef "vect"];[cFunDef "main"; cFor "i"]]];

  (* Testing of marks *)
  (* TODO: test
  let mark = Ast.next_mark() in
  Generic.add_mark mark [ nbMulti; cFor "i" ];
  Transfo.iter_targets [cMark mark_bindings] (fun (p:path) ->
     Generic.add_mark (Ast.next_mark()) (target_of_path p)
    );
  Generic.remove_marks [ cMark mark ];
  Generic.remove_marks [ cMarkAll  ];
  *)
  (* TODO: print the marks when printing the code, using a similar code as for annotation
      like /*<@*/ for show,  would be for example   /*@43*/ *)

(* TODO:
let Generic_core.add_mark (mark : mark) (t : trm) : trm =
  let left_decoration = "/*[" ^ string_of_int mark ^ "]*/" in
  let right_decoration = "" in
  {t with annot = t.annot @ [Highlight (left_decoration, right_decoration)]}


  *)


    (* Testing the [target_of_path] function *)
  !! let paths = ref [] in
  Target.apply_on_targets (fun t p -> paths := p :: !paths; t) [ nbMulti; cFor "i" ];
  let targets = List.map target_of_path !paths in
  let tg1,tg2 =
    match targets with
    | target1:: target2 :: [] -> target1,target2
    | _ ->  assert false
    in
  show tg1;
  show tg2;
)


