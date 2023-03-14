include Marks_basic

open Target
open Ast

let with_fresh_mark (f : mark -> unit) : unit =
  let m = Mark.next () in
  f m;
  remove m [nbAny; cMark m]

let with_fresh_mark_on (p : path) (f : mark -> unit) : unit =
  with_fresh_mark (fun m ->
    add m (target_of_path p);
    f(m)
  )