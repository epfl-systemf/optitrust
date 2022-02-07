open Optitrust
open Target

(* let _ = Run.doc_script_cpp (fun _ ->
  !! Instr_basic.read_last_write ~write:[sInstr "a ="] [cVarDef "b"; dBody];
  )
"
int main() {
  int a = 5;
  a = 6;
  int b = a;
}
" *)

let _ = Run.script_cpp (fun _->

    !! Instr_basic.read_last_write ~write:[cWrite ~rhs:[cInt 7] ()] [cRead ~addr:[cVar "x" ] ()];
    (* LATER: fix printing of int a ;  should have no space *)
    (* !! Instr.view_subterms ~constr:(sInstr "int a ;") [dRoot];*)
    show [tBefore;cVarDef "a"];
    show [ tAfter; (* sInstr "int a ;" *) cVarDef "a"] ;
    (*
    !! Instr_basic.read_last_write ~write:[sInstr "t[0] ="] [sInstr "= t[0]"; dRHS];*)
)
