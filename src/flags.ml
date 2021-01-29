(* debug printing*)
let verbose : bool ref = ref false

(* interprete dump as dump_trace in transformation scripts *)
let full_dump : bool ref = ref false

(*
  force printing/parsing at each transformation step
 *)
let repeat_io : bool ref = ref false

let spec =
  Arg.align [
     ("-verbose", Arg.Set verbose, " activates debug printing");
     ("-dump-trace", Arg.Set full_dump, " dump ouptputs the full trace in " ^
                                          "transformation scripts");
     ("-repeat-io", Arg.Set repeat_io, " print/parse the current program at " ^
                                         "each transformation step")
    ]
