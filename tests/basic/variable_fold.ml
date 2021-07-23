open Optitrust
open Target

let _ = Run.script_cpp( fun _ ->
  !! Variable.fold [cVarDef "s1" ];
  !! Variable.fold [cVarDef "s2" ];
  !! Variable.fold [cVarDef "a" ];
  (* TODO: demo with ~at: *)
  (* LATER:
     at the basic level,
        Variable.fold_basic  works for both const and non-const
     at a high level:
     - Variable.fold  works for const, fails for nonconst
     - Variable.fold ~nonconst:true    : to force it working for nonconst
   *)
)