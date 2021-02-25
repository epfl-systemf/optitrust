open ScriptTools

let _ = 
    run(
        fun _ -> 
        set_init_source "coloring.cpp";
        loop_transform [cFor ~init:[cVarDef ~name:"i" ()] ()] "C" ;
        loop_transform [cFor ~init:[cVarDef ~name:"j" ()] ()] "C" ;
        loop_tile [cFor ~init:[cVarDef ~name:"x" ()] ()] "B";
        dump()
    )