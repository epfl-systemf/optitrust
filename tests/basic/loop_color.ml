open Optitrust
open Target

let _ = Run.script_cpp (fun _ ->
  !! Loop.color "C" ~index:"ci" [cFor "i"];
  !! Loop.color "C" [cFor "j"];
  (* LATER: in the future, we will check that "ci" is not conflicting
     with an existing variable *)
)
