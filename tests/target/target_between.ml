open Optitrust
open Target


let _ = Run.script_cpp (fun () ->
  (** There should be exactly one result to each of the commands;
      if it is not the case, we'll get an error. *)
  (* Before *)
  show [cIf (); dThen];
  show [ tBefore; cVarDef "r1" ];
  show [ tBefore; cVarDef "r2" ];
  show [ tBefore; cVarDef "m1" ];
  show [ tBefore; cVarDef "m2" ];

  (* After *)
  show [ tAfter; cVarDef "r1" ];
  show [ tAfter; cVarDef "r2" ];
  show [ tAfter; cVarDef "m1" ];
  show [ tAfter; cVarDef "m2" ];

  (* First *)
  show [ tFirst; cFor "i"; cStrict; dBody ];
  show [ tFirst; dElse ];

  (* Last *)
  show [ tLast; dElse ];
  show [ tLast; cFor "i"; cStrict; dBody];
  show [ tLast; dThen ];
  show [ tLast; dElse ];

  (* Nested paths *)
  show [ tLast; cFor_c"i"; dBody; dThen ];
  (* Top level paths *)
  show [tBefore; cTopFunDef "main"];
  show [tAfter; cTopFunDef "main"];

)


