open Optitrust


let _ = run_unit_test (fun () ->
  let show = show_target in
  set_repeat_io false;
  show [cMulti; cVar "u"];
  (* show [ cTypDef "intstar" ]; *)
  (* show [ cMulti; cFunDef "f" ]; *)
  (*  show [ cTopFun "f" ];*)
  (*    show [ cTopFun "main" ];*)
)
