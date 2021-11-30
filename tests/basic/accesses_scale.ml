open Optitrust
open Target
open Ast

let _ = Run.doc_script_cpp (fun _ ->
    !! Accesses_basic.scale (trm_double 5.0) [cRead ~addr:[cVar "x"]()]
  )
"
int main() {
  int x = 2;
  int y = x;
}
"

let _ = Run.script_cpp (fun _ ->

  !! Accesses_basic.scale (trm_double 5.0) [cOr [[cCellWrite ~base:[cVar "t"] ~index:[cVar "i"] ()];[cCellRead ~base:[cVar "t"] ~index:[cVar "i"] ()]]]
)
