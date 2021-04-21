open Optitrust


let _ = 
  run
    ( fun () ->
      set_init_source "pic_demo.cpp";
      (*
      inline_decl ~delete_decl:false ~decl_path:[>cTopFun ~name:"v_add" ()]~fun_result:"myres" (*~keep_labels:true *) ();
      inline_decl ~delete_decl:true ~decl_path:[cVarDef ~name:"myres"()] ();
      inline_decl ~delete_decl:true ~decl_path:[cVarDef ~name:"v1_0"()] ();
      inline_decl ~delete_decl:true ~decl_path:[cVarDef ~name:"v2_0"()] ();*)
      detach_expression [cVarDef ~name:"speed2"()] ~keep_label:false; 
            (* TODO: later: detach all of a given type *)
      inline_decl ~delete_decl:false ~decl_path:[cTopFun ~name:"vect_add" ()] ~inline_at:[[cTopFun ~name:"main"()]]();
      (*detach_expression [cVarDef ~name:"pos2"()] ~keep_label:false; 
      inline_decl ~delete_decl:false ~decl_path:[cTopFun ~name:"vect_add" ()] (*~keep_labels:true *) (); *)
      inline_decl ~delete_decl:true ~decl_path:[cVarDef ~name:"res" ()] ~inline_at:[[cTopFun ~name:"main"()]]();
      make_explicit_record_assignment ~struct_name:"vect" [cSet ~lhs:[cVar ~name:"speed2"()]  ()]; 
      detach_expression [cVarDef ~name:"v2"()] ~keep_label:false; 
      inline_decl ~delete_decl:false ~decl_path:[cTopFun ~name:"vect_mul" ()] ~inline_at:[[cTopFun ~name:"main"()]]();
      inline_decl ~delete_decl:true ~decl_path:[cVarDef ~name:"res" ()] ~inline_at:[[cTopFun ~name:"main"()]]();
      inline_decl ~delete_decl:true ~decl_path:[cVarDef ~name:"v1"()] ~inline_at:[[cTopFun ~name:"main"()]]();
      inline_decl ~delete_decl:true ~decl_path:[cVarDef ~name:"v"() ] ~inline_at:[[cTopFun ~name:"main"()]] ();
      inline_decl ~delete_decl:true ~decl_path:[cVarDef ~name:"d"()]  ~inline_at:[[cTopFun ~name:"main"()]] ();
      inline_record_access ~field:"x" ~var:"v2" ();     
      inline_record_access ~field:"y" ~var:"v2" ();     
      inline_record_access ~field:"z" ~var:"v2" ();
      (* delete_decl  *)
      detach_expression [cVarDef ~name:"pos2"()] ~keep_label:false; 
      inline_decl ~delete_decl:true ~decl_path:[cTopFun ~name:"vect_add" ()] (*~keep_labels:true *) ();
      inline_decl ~delete_decl:true ~decl_path:[cVarDef ~name:"res" ()] (*~keep_labels:true *) ();
      make_explicit_record_assignment ~struct_name:"vect" [cSet ~lhs:[cVar ~name:"pos2"()]  ()]; 
      
      inline_struct ~struct_name:"particle" ~struct_fields:["speed"] ();
      



      
      (* show_ast [cVarDef ~name:"v2"()] ; *)
       (*show_path[cFun ~name:"main" (); cVarDef ~name:"v1" ()] ;
        TOOD: need to fix paths
      inline_decl ~delete_decl:true ~decl_path:[cFun ~name:"main" (); cVarDef ~name:"v1" ()] (*~keep_labels:true *) ();
      inline_decl ~delete_decl:true ~decl_path:[cVarDef ~name:"d" ()] (*~keep_labels:true *) ();
      inline_decl ~delete_decl:true ~decl_path:[cVarDef ~name:"v" ()] (*~keep_labels:true *) ();
      show_path [cFor ~init:[cVar ~name:"idParticle" ()] ()]; *)
      
       (* show_ast [cVarDef ~name:"v2"()] ; *)
      (* show_path [cSet ~lhs:[cVar ~name:"u"()] ()] ~debug_ast:true; *)
      (* show_ast [cVarDef ~name:"v1"()] ; *)
      (* show_ast [cSet ~lhs:[cVar ~name:"speed2"()]()];  *)
      
      
      
      
      
      
      (* show_path [cApp ~args:[cVar ~strict:true ~name:"p2" ()] ~validate:(List.mem true) ()] ~debug_ast:true; *)
      inline_decl ~delete_decl:false~decl_path:[cFun ~name:"bag_push" ()] ~inline_at:[[cTopFun ~name:"main"()]]();
      
      (* TODO: Fix the error with variable bound *)
      
      (* make_explicit_record_assignment ~struct_name:"particle" [cApp ~args:[cVar ~strict:true ~name:"p_0" ()] ~validate:(List.mem true) ()]; *)
      (* TODO: if struct_name is not provided, it is infered from the type of the LHS or RHS of the assignement *)
      (* LATER: also genertate a show_path on the path involved in the last transformation. *)
(* make_explicit_record_assignment *)
      dump()
    )


    (*
    TODO

     inline_seq_core : takes a path to the outer sequence, and the index of the statement corresponding to the inner sequence
     inline_seq:  takes a path to the inner sequence
       and computes the directions to it, it should end with "path_to_outerseq@[Dir i]" (use list_extract_last , and general version is [list_extract_lastn n l]), then call inline_seq_core.

     call_inline ~decl:path ~name:string (interpreted as toplevel function) ~path:path (to specific call sites) ~keep_decl ~arg_names:(string list) ~res_name:(..) ~label:string=""
     const_inline ~decl:path ~name:string (interpreted as cVarDef) ~path:path (to specific occurences)


    const int a = 3


    a // replacing by 3 is always ok


1a.    
      int a = 3 // if not a constant
      b = a 
1b
      -> b = 3    (inline_decl "a")

1a is equivalent to 2a: // take this view when not a constant
      int a;
      a = 3;
      b = a
->    
      int a;
      a = 3;
      b = 3

  transformation suggested:  read_inline_last_write

     set (ptr, v)    // path1

     ...

     get (ptr)      //path2


     idea is to replace the get(ptr) by v
     the user takes responsibility for the fact that there is no other write in the middle


------------------------

   let inlining_core ...
     // current code, does no simplification


   let inlining ?(simpl:bool=true) =
     inline_core ~labels:["block";"arg";"res"].
     let ts = get the terms at label "arg_*"
     for each result you have a line of "arg_i: int a2 = v3;"
        if v3 looks like a constant or a variable,
        invoke inline_decl at path  "arg_i", strictly below
      



   
    *)

          (* for x = 0; x < X;
          for y = 0; y < Y
        flattening: assumes rectangular area

          for k = 0; k<X*Y;k++ 
             x = k / Y
             y = k % Y
       *)
