open Optitrust
open Target

(* TODO
let _ = Run.doc_script_cpp (fun _ ->
  !! Rewrite.equiv_at "int x; ==> 3 + x * 4 == 4 * x + 3;" [cVarDef "b"; dBody];
  )
"
int main() {
  int a = 2;
  int b = 3 + a * 4;
}
"
*)

let _ = Run.script_cpp (fun _ ->

  !! Rewrite.equiv_at "double a, b; int k; ==> a + k * b == b * k  + a;" [cWriteVar "res"; dRHS];
  !! Rewrite.equiv_at "double a; int k; ==> a + k * a == (k + 1) * a;" [cWriteVar "res1"; dRHS];
  !! Rewrite.equiv_at "double a; int k; ==> a + k * a == (k + 1) * a;" [cVarDef "res2"; dBody];
)
