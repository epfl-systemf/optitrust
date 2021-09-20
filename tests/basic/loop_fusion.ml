open Optitrust
open Target

let _ = Run.script_cpp ( fun _ ->
  
  !! Sequence_basic.intro ~label:"tofusion" 3 [cFunDef "fusion_on_block"; cFor "i" ~body:[sInstr "t[i]"]];
  !! Loop_basic.fusion_on_block [cLabel "tofusion"];  
)
