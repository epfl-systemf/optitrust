open Optitrust
open Target

let _ = Run.script_cpp (fun () ->
  !! Loop_basic.hoist "x_step" [cVarDef "x"];
  !! Loop_basic.hoist "z_step" [cVarDef "z"];
)

(* LATER: in combi level, add
      Loop.hoist = detach if needed + Loop.hoist_without_detach

    if before
      int x = 2
    first need to detach it.
      int x
      x = 2
    by detaching first of all.
*)
