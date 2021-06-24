open Optitrust
open Target

let _ = Run.script_cpp ( fun _ ->
    !! Label.add "start" [cVarDef "x"] ;
    !! Label.add "loop" [cFor "i"];
    !! Label.add "cond" [cIf ~then_:[cVar "x"] ()] ;
    !! Label.add "incr_1" [cIf (); cVar "x"];
    !! Label.add "incr_2" [cIf (); sInstr "i++" ];    
    !! Label.add "stop" [sInstr "return"];
)
