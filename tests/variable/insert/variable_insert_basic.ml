open Optitrust
open Prelude

let _ = Run.script_cpp (fun _ ->

  !! Variable_basic.insert ~const:true ~name:"a" ~typ:(ty "int") ~value:(lit "300") [ tAfter; cTypDef "vect"];
  !! Variable_basic.insert ~reparse:true ~name:"b" ~typ:(ty "int") ~value:(lit "500") [ tAfter; cTypDef "vect"];
  !! Variable_basic.insert ~name:"c" ~typ:(ty "int") [tAfter ; cTypDef "vect"];

)
