open Optitrust
open Run

(* TODO: Noy yet implemented *)

let _ = 
    run 
    ( fun _ -> 
        set_init_source"make_implicit_record_assignment.cpp";

        make_explicit_record_assignment [cVarDef "b"] ~struct_name:"vect";

        make_implicit_record_assignment [cVarDef "b"] ~struct_name:"vect" ; 
        
        make_implicit_record_assignment [cVarDef "d"] ~struct_name:"vect" ; 
        
        Generic.var_init_detach [cVarDef "b"] ~keep_label:true; 
        make_explicit_record_assignment [cLabel "detached";cBody()] ~struct_name:"vect";
        
        make_explicit_record_assignment [cVarDef "b"] ~struct_name:"vect";
        
        make_explicit_record_assignment [cVarDef "b"] ~struct_name:"vect";
        delete_label "detached";

        dump()
    )
