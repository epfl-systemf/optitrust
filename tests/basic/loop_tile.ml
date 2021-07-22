open Optitrust
open Target

let _ = Run.script_cpp (fun _ ->
  !! Loop.tile "2" ~index:"bx" ~bound:TileBoundDivides [cFor "x"];
  !! Loop.tile "2" ~bound:TileBoundMin [cFor "y"];
  !! Loop.tile "2" ~bound:TileBoundAnd [cFor "z"];

  !! Loop.tile "2" ~index:"ba" ~bound:TileBoundDivides [cFor "a"];
  !! Loop.tile "2" ~bound:TileBoundMin [cFor "b"];
  !! Loop.tile "2" ~bound:TileBoundAnd [cFor "c"];
  (* TODO: fix *)


  (* TODO: rename to i,j,k *)


  (*
  type tile_bound = TileBoundMin | TileBoundAnd | TileBoundDivides

     Loop.tile ~bound:TileBoundDivides


  TODO: MOVE THESE TO THE ML FILE DOCUMENTATION

  for (int x = 0; x < 10; x++) {
    -->
  for (int bx = 0; (bx < 10); bx += 2) {
    for (int x = bx; (x < bx + 2); x++) {


     Loop.tile  ~divides:TileBoundMin // make the default

  for (int i = 0; i < 9; i++) {
    -->
  for (int bi = 0; (bi < 9); bi += 2) {
    for (int i = bx; i < min(bi + 2, 9); i++)


     Loop.tile  ~divides:TileBoundAnd

  for (int bi = 0; (bi < 9); bi += 2) {
    for (int i = bx; i < bi + 2 && i < 9; i++)



  --- With step, the parameter of Loop.tile tells how many
  iterations of the body should be performed (here we assume 2)

  for (int y = 0; y < 80; y += 3) {


  for (int by = 0; (by < 100); by += 2 * 3) {
    for (int y = by; ...like before based on (y < by + 2 * 3) ...; y += 3) {
  *)
)
