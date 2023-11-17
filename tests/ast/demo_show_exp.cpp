#include "../../include/optitrust.h"

int main() {
  int a, b;
  int x = 3;
  x--;
  for (int i = 0; i < 3; i++) {
    __sequentially_modifies("x ~> Cell");
    /*@mymark2, mymark1*/ x /*mymark2, mymark1@*/++;
  }
}
/*
CAPTURED STDOUT:
AST1: 

  int main (  )  {
  int a, b;
  int x = 3;
  x--;
  for ( int i = 0; i < 3; i++ ) {
    __sequentially_modifies("x ~> Cell");
    x++;
  }
}
AST2: 

  int main (  )  {
  int a, b;
  int x = 3;
  x--;
  for ( int i = 0; i < 3; i++ ) {
    __sequentially_modifies("x ~> Cell");
    x++;
  }
}
var-trm-internal:   int* x = new int( 3 )
for-trm-internal: for ( int i = 0; i < 3; i++ ) { x++; }
for-trm-internal-desc: Trm_for ( i#351,
  Trm_val ( Val_lit ( Lit_int 0 ) ),
  Up,
  Trm_val ( Val_lit ( Lit_int 3 ) ),
  Post_inc,
  false,
  Trm_seq [
    Trm_apps ( Trm_val ( Val_prim ( Prim_unop Unop_post_inc ) ),
      [ Trm_var( Var_mutable,x#350) ] )
  ] )
desc: Trm_for
typ: target resolves to 2 paths
[occ #1] int
[occ #2] int
marks: [mymark2; mymark1]
cstyle-item: [ Stackvar ]
annot: {trm_annot_attributes = [  ]; trm_annot_marks = [  ];
  trm_annot_stringrepr = None; trm_annot_pragma = [  ]; trm_annot_cstyle = [

]; trm_annot_files = [ Main_file ]}

*/
