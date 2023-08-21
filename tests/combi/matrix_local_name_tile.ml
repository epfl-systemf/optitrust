open Optitrust
open Target
open Syntax

let _ = Run.doc_script_cpp (fun _ ->
  let a = find_var_in_current_ast "a" in
  !! Matrix.local_name_tile a ~into:"b" [(trm_int 3, trm_int 4)] [cFor "i"];
)

"
#include \"../../include/optitrust.h\"

int main (){
  int* a = (int*) CALLOC1(10, sizeof(int));
  for (int i = 3; i < 7; i++) {
    a[MINDEX1(10, i)];
  }
  MFREE1(10, a);
}
"

let _ = Run.script_cpp (fun _ ->
  let fv = find_var_in_current_ast in

  let tile offset size = (trm_int offset, trm_int size) in
  !! Matrix.local_name_tile (fv "a") [tile 0 10; tile 2 8; tile 0 4] [cFor ~body:[cArrayWrite "a"] "i"];
  !! Matrix.local_name_tile ~delete:true (fv "b") ~into:"y" [tile 0 10; tile 2 8; tile 0 4] ~alloc_instr:[cVarDef "b"] [cFor ~body:[cArrayWrite "b"] "j"];
  !! Matrix.local_name_tile ~delete:true (fv "c") ~into:"z" [(var "i", trm_int 1); tile 2 8; tile 0 4] ~alloc_instr:[cVarDef "c"] [cFor ~body:[cArrayWrite "c"] "j"];
)
