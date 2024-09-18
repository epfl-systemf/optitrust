open Optitrust
open Prelude
open Option

let _ = Flags.pretty_matrix_notation := false
let _ = Flags.disable_resource_typing ()

let int = trm_int

let _ = Run.script_cpp (fun () ->
  !! (fun () -> ())();
  (*!! Arrays_basic.tile (find_var_in_current_ast "B") [cVarDef "sum"];*)
  (*!! Nobrace_transfo.remove_after ~remove:true (fun _ -> Variable.bind_multi ~dest:[cIf (); tBefore] "ofs" [sExpr "tile_j * 4 + j + r"]);
  !! Matrix_basic.elim_mindex [cMindex ~args:[[cVar "IC"]; [cVar "IW"]; [cVar "k"]; [cVar "ofs"]] ()];
  !! Matrix_basic.elim_mindex [nbMulti; cMindex ~args:[[cInt 4]; [cInt 4]; [cVar "j"]; [cVar "r"]] ()];
  !! Variable.bind_multi ~dest:[cFor "r"; tBefore; occFirst] "dtile_ofs" [sExpr "0 + j * 4 + r"];*)
)
