(* option: use option to add open at compile time *)
open ScriptTools

let _ =
  run
    (fun () ->

      set_init_source "test_split.cpp";
      add_label "return_instr" [cInstrSubstr "return"];

      split_loop ~keep_labels:false [cInstrSubstr ~regexp:true "^x ="];
      
      dump ()

    )
