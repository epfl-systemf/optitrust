open Optitrust
open Target

let _ = Run.script_cpp (fun _ ->

   show [cVarDef "p"];
   !! Variable.inline [cVarDef "p"];
)
