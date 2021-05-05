open Optitrust

let _ = run_unit_test (fun () ->
  (** There should be exactly one result to each of the commands;
      if it is not the case, we'll get an error. *)
  let show = Tr.target_show in

  (* Types *)
  show [ cTypDef "vect" ];
  show [ cTypDef "intstar" ];

  (* Constants *)
  show [ cInt "8"; ];

  (* Var/fun occurences *)
  show [ cVar "u" ];
  show [ cVar "r2" ];
  show [ cVar "f" ];
  show [ cVar "g" ];

  (* Loops *)
  show [ cFor "i" ];
  show [ cFor "j" ];
  show [ cFor ~cond:[cStr "j < 5"] "" ];

  (* Abort *)
  show [ cBreak ];
  show [ cContinue ];
  show [ cReturn ];

  (* Labels *)
  show [ cLabel "lbl1" ];
  show [ cLabel "lbl2" ];

  (* Calls *)
  show [ cCall "f" ];
  show [ cCall ~args:[[cInt 2]] "" ];

  (* Var/Fun definitions *)
  show [ cDef "f" ];
  show [ cDef "s" ];
  show [ cDef "p2" ];
  show [ cFun "main" ];
  show [ cFun "f" ];
  show [ cFun ~args:[[cTrue]; [cDef "varg"]] "" ];
  show [ cFun ~argspred:((fun i -> [cTrue]),(fun bs -> List.length bs = 2)) "" ];

)



(* LATER: smart constructors for checking calls to builtin operations such as get/set/compare/incr, etc *)

(* LATER: show [ cFunDef ~args:[[cTrue]; [cOfTyp "vect*"]] "" ]; *)

(* LATER: match typedef using a function over the body of the type definition *)
(* LATER: match a typedef struct using of a function over the list fields [(var*typ)list->bool] *)

(* LATER: match types using a function of their list of fields *)

(* LATER: Str/Regexp -- i still not to work on the specification of what is expected to match or not match
  show [ cStr "+= 2" ];
  show [ cStr "j <" ];
  show [ cStr "vect v2" ];
  show [ cStrFull "int r = 3;" ]; (* with or without the ; ? *)
  show [ cStrFull "i++" ];
  show [ cRegexp "+=" ];
  show [ cRegexp ~sub:false "+=" ];
  show [ cRegexp "f\\(.\\)" ];
  *)



