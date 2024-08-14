open Optitrust
open Prelude
open Option

let _ = Flags.pretty_matrix_notation := false
let _ = Flags.disable_resource_typing ()

let _ = Run.script_cpp (fun () ->
  !! (fun () -> ())()

)
