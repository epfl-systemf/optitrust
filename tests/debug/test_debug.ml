open Optitrust
open Target

let _ = Flags.dump_ast_details := true

let _ = Run.script_cpp (fun _ ->

   !! Variable_basic.ref_to_pointer [cVarDef "x"];
   
)

