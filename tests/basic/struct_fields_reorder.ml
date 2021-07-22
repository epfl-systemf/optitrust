
open Optitrust
open Target

let _ = Run.script_cpp (fun _ ->
    !! Struct.fields_reorder ~move_before:"x" ["m";"z"] [cTypDef "obj"];
    !! Struct.fields_reorder ~move_after:"y" ["z"] [cTypDef "obj"];      
    !! Struct.fields_reorder ~move_after:"z" ["m"] [cTypDef "obj"];      
)
