open Optitrust
open Target

let _ = Run.script_cpp (fun _ ->

  !! Variable.reuse ~reparse:true "x" [cVarDef "y"];
)
