open Optitrust
open Target
open Ast

let _ =  Run.script_cpp ( fun () ->
  !! Variable.delocalize ~loop_index:"k" ~mark:"section_of_interest" ~old_var:"a" ~new_var:"x" ~var_type:(typ_constr "T") ~array_size:"N" ~dl_ops:(Delocalize_arith (Lit_int 0, Binop_add) ) [cFor "i"];
  !! Trace.alternative (fun () ->
    !! Variable.local_other_name  ~var_type:(typ_constr "T") ~old_var:"a" ~new_var:"x" [cFor "i"];
    !! Variable_basic.delocalize  ~array_size:"N" ~dl_ops:(Delocalize_arith (Lit_int 0, Binop_add)) [cMark "section_of_interest"] ;
    !! ());
  !! Trace.alternative (fun () ->
    !! Variable.local_other_name  ~var_type:(typ_constr "T") ~old_var:"a" ~new_var:"x" [cFor "i"];
    !! Variable_basic.delocalize  ~array_size:"N" ~dl_ops:(Delocalize_obj ("clean", "transfer")) [cMark "section_of_interest"] ;
    !! ());
)
