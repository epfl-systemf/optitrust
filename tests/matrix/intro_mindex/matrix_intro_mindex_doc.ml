open Optitrust
open Target
open Prelude


let _ = Run.script_cpp (fun () ->

  !! Matrix.intro_mindex (var "N") [cVarDef "p"];

)
