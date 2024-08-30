open Optitrust
open Prelude
open Option

let _ = Flags.pretty_matrix_notation := true
let _ = Flags.disable_resource_typing ()

let _ = Run.script_cpp (fun () ->
  !! Loop.reorder ~order:["k"; "j"] [cFor "j"];
)
