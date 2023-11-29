open Ast
open Trm
open Typ
open Mark

(**
  A [formula] is a [trm] that corresponds to a logical formula. Evaluating inside a formula cannot have side effects and always terminates.
  Formulas may refer to pure variables, bound by a [resource_set] or by the program (program variables are constant due to our AST encoding).
*)

(** An optionally named resource item produced by parser.
    TODO: remove. *)
type contract_resource_item = var option * formula

let new_hyp = new_var

(** Number used to generate variable names for resources. *)
let next_hyp_id = Tools.fresh_generator ()

(** Returns a variable with a generated name.

  TODO: separate pure ($) and linear (#). *)
let new_anon_hyp (): hyp =
  let hid = next_hyp_id () in
  new_hyp (sprintf "#%d" hid)

(* TODO: should be new_var_like and maybe useful elsewhere? *)
let new_hyp_like (h: hyp): hyp =
  new_hyp ~qualifier:h.qualifier h.name

(** _HasModel(p, Cell) <=> p ~> Cell *)
let var_has_model = toplevel_var "_HasModel"
let trm_has_model = trm_var var_has_model

(** Primitive function that constructs a read only resource. *)
let var_read_only = toplevel_var "_RO"
let trm_read_only = trm_var var_read_only

(** Primitive function that constructs an uninit resource. *)
let var_uninit = toplevel_var "_Uninit"
let trm_uninit = trm_var var_uninit

(** Primitive type of fractions. *)
let var_frac = toplevel_var "_Fraction"
let trm_frac = trm_var var_frac

(** Creates a new fraction variable.
    It ranges on fraction values from \]0; 1\]. *)
let new_frac (): var * resource_item =
  let frac_hyp = new_anon_hyp () in
  (frac_hyp, (frac_hyp, trm_frac))

(** The fraction representing having it all. *)
let full_frac = trm_int 1

(** All formulas should have this annotation. *)
let formula_annot = {trm_annot_default with trm_annot_cstyle = [ResourceFormula]}

(** Tries to embed a program term within formulas.
    Pure and total terms can be successfully embedded, according to built-in whitelist. *)
let rec formula_of_trm (t: trm): formula option =
  (* LATER: Extensible list of applications that can be translated into formula.
     OR A term can be embedded in a formula if it is pure (no linear resources).
     Does [t] change or is this just checking whether trm_is_formula?
     I.e. are the trm and formula languages intersecting or separate?
     *)
  let open Xoption.OptionMonad in
  match t.desc with
  | Trm_val _ | Trm_var _ -> Some t
  | Trm_apps (fn, args, _) ->
    let* f_args = try Some (List.map (fun arg -> Option.get (formula_of_trm arg)) args) with Invalid_argument _ -> None in
    begin match trm_prim_inv fn with
      | Some Prim_binop Binop_add
      | Some Prim_binop Binop_sub
      | Some Prim_binop Binop_mul
      | Some Prim_binop Binop_div (* TODO: think hard about 'div' totality *)
      | Some Prim_binop Binop_mod
      | Some Prim_binop Binop_array_access
          -> Some (trm_apps fn f_args)
      | Some _ -> None
      | None ->
        begin match trm_var_inv fn with
          | Some fv when Matrix_trm.mindex_var_inv fv <> None
              -> Some (trm_apps fn f_args)
          | _ -> None
        end
    end
  | _ -> None

(* -------- SMART FORMULA CONSTRUCTORS, INVERTERS and COMBINATORS -------- *)

let formula_fun =
  trm_fun ~annot:formula_annot

let formula_model (x: trm) (model: formula): formula =
  trm_apps ~annot:formula_annot trm_has_model [x; model]

let formula_model_inv (t: formula): (trm * formula) option =
  match trm_apps_inv t with
  | Some (fn, [tx; tf]) ->
    begin match trm_var_inv fn with
    | Some fnv when var_eq var_has_model fnv -> Some (tx, tf)
    | _ -> None
    end
  | _ -> None

let formula_read_only ~(frac: formula) (res: formula) =
  trm_apps ~annot:formula_annot trm_read_only [frac; res]

let formula_uninit (inner_formula: formula): formula =
  trm_apps ~annot:formula_annot trm_uninit [inner_formula]

let var_cell = toplevel_var "Cell"
let trm_cell = trm_var var_cell

let formula_cell (x: var): formula =
  formula_model (trm_var x) trm_cell

let var_group = toplevel_var "Group"
let trm_group = trm_var var_group
let var_range = toplevel_var "range"
let trm_range = trm_var var_range

let formula_matrix (m: trm) (dims: trm list) : formula =
  let indices = List.mapi (fun i _ -> new_var (sprintf "i%d" (i+1))) dims in
  let inner_trm = formula_model (Matrix_trm.access m dims (List.map trm_var indices)) trm_cell in
  List.fold_right2 (fun idx dim formula ->
    trm_apps ~annot:formula_annot trm_group [trm_apps trm_range [trm_int 0; dim; trm_int 1]; formula_fun [idx, typ_int ()] None formula])
    indices dims inner_trm

module Pattern = struct
  include Pattern

  let formula_model f_var f_model =
    trm_apps_specific_var var_has_model (f_var ^:: f_model ^:: nil)

  let formula_read_only f_frac f_formula =
    trm_apps2 (trm_var (var_eq var_read_only)) f_frac f_formula

  let formula_uninit f_formula =
    trm_apps1 (trm_var (var_eq var_uninit)) f_formula

  let formula_group f_range f_group_body =
    trm_apps_specific_var var_group (f_range ^:: f_group_body ^:: nil)

  let formula_range (f_begin: 'a -> trm -> 'b) (f_end: 'b -> trm -> 'c) (f_step: 'c -> trm -> 'd) =
    trm_apps_specific_var var_range (f_begin ^:: f_end ^:: f_step ^:: nil)
end

type read_only_formula = { frac: formula; formula: formula }
let formula_read_only_inv (formula : formula): read_only_formula option =
  Pattern.pattern_match_opt formula [
    Pattern.(formula_read_only !__ !__) (fun frac formula -> { frac; formula })
  ]

(** Applies a function below a read only wrapper if there is one,
    otherwise simply applies the function. *)
let formula_map_under_read_only (f_map: formula -> formula) (formula: formula) =
  match formula_read_only_inv formula with
  | Some { frac; formula } ->
    formula_read_only ~frac (f_map formula)
  | None -> f_map formula

let formula_uninit_inv (formula: formula): formula option =
  Pattern.pattern_match_opt formula [
    Pattern.(formula_uninit !__) (fun f -> f);
  ]


(** Applies a function below an uninit wrapper if there is one,
    otherwise simply applies the function. *)
let formula_map_under_uninit (f_map: formula -> formula) (formula: formula) =
  match formula_uninit_inv formula with
  | Some formula -> formula_uninit (f_map formula)
  | None -> f_map formula

(** Applies a function below a resource mode wrapper if there is one,
    otherwise simply applies the function.

  Current resource mode wrappers include read only and uninit, they float to the top of formulas.
 *)
let formula_map_under_mode (f_map: formula -> formula): formula -> formula =
  formula_map_under_read_only (formula_map_under_uninit f_map)

let formula_loop_range ((_, tfrom, dir, tto, step, _): loop_range): formula =
  if dir <> DirUp then failwith "formula_loop_range only supports DirUp";
  trm_apps trm_range [tfrom; tto; loop_step_to_trm step]

let formula_group_range ((idx, _, _, _, _, _) as range: loop_range) =
  formula_map_under_mode (fun fi ->
    let range_var = new_var ~qualifier:idx.qualifier idx.name in
    let fi = trm_subst_var idx (trm_var range_var) fi in
    trm_apps ~annot:formula_annot trm_group [formula_loop_range range; formula_fun [range_var, typ_int ()] None fi]
  )

let formula_matrix_inv (f: formula): (trm * trm list) option =
  let open Xoption.OptionMonad in
  let rec nested_group_inv (f: formula): (formula * var list * trm list) =
    Pattern.pattern_match f [
      Pattern.(formula_group (formula_range (trm_int (eq 0)) !__ (trm_int (eq 1))) (trm_fun (pair !__ __ ^:: nil) !__))
        (fun dim idx inner_formula ->
          let inner_formula, indices, dims = nested_group_inv inner_formula in
          (inner_formula, idx::indices, dim::dims)
        );
      Pattern.__ (f, [], [])
    ]
  in
  let inner_formula, indices, dims = nested_group_inv f in

  let* location, cell = formula_model_inv inner_formula in
  let* cell_candidate = trm_var_inv cell in
  let* () = if var_eq cell_candidate var_cell then Some () else None in
  let* matrix, mindex_dims, mindex_indices = Matrix_trm.access_inv location in
  let* () = if List.length mindex_dims = List.length dims then Some () else None in
  let* () = if List.for_all2 are_same_trm mindex_dims dims then Some () else None in
  if List.for_all2 (fun mindex_idx idx ->
      match trm_var_inv mindex_idx with
      | Some idx_arg when var_eq idx idx_arg -> true
      | _ -> false
    ) mindex_indices indices
  then Some (matrix, dims)
  else None

let var_fun_type = toplevel_var "_Fun"

let formula_fun_type (targ: trm) (tres: trm) =
  trm_apps (trm_var var_fun_type) [targ; tres]

let formula_assert_eq = toplevel_var "__assert_eq"
let formula_assert_neq = toplevel_var "__assert_neq"
let formula_assert_lt = toplevel_var "__assert_lt"
let formula_assert_gt = toplevel_var "__assert_gt"
let formula_assert_leq = toplevel_var "__assert_leq"
let formula_assert_geq = toplevel_var "__assert_geq"

let formula_cmp (cmp: var) (a: formula) (b: formula): formula =
  trm_apps ~annot:formula_annot (trm_var cmp) [a; b]

let formula_eq = formula_cmp formula_assert_eq
let formula_neq = formula_cmp formula_assert_neq
let formula_lt = formula_cmp formula_assert_lt
let formula_gt = formula_cmp formula_assert_gt
let formula_leq = formula_cmp formula_assert_leq
let formula_geq = formula_cmp formula_assert_geq

let var_checked = toplevel_var "checked"
let formula_checked = trm_var var_checked

(* ---- RESOURCE-RELATED TRM CONSTRUCTORS and INVERTERS ---- *)
(* TODO: Resource_trm ? *)

(** Primitive function that turns off resource computation and checking from that point onwards.
    Resource computation resumes upon meeting an outer postcondition which will be admitted. *)
let __admitted = toplevel_var "__admitted"

let trm_ghost_inv t =
  Pattern.pattern_match t [
    Pattern.(trm_apps !__ nil !__) (fun ghost_fn ghost_args ->
      if not (trm_has_cstyle GhostCall t) then raise Pattern.Next;
      Some { ghost_fn; ghost_args }
    );
    Pattern.(__) None
  ]

let ghost_begin = toplevel_var "__ghost_begin"
let ghost_end = toplevel_var "__ghost_end"
let typ_ghost_fn: typ = typ_const (typ_constr ([], "__ghost_fn") ~tid:(-1))
let with_reverse = toplevel_var "__with_reverse"

let ghost_pair_fresh_id = Tools.fresh_generator ()
let generate_ghost_pair_var ?name () =
  match name with
  | Some name -> new_var name
  | None -> new_var (sprintf "__ghost_pair_%d" (ghost_pair_fresh_id ()))

let trm_ghost_begin (ghost_pair_var: var) (ghost_call: ghost_call) : trm =
  trm_add_cstyle GhostCall (trm_let Var_immutable (ghost_pair_var, typ_ghost_fn)
  (trm_apps (trm_var ghost_begin) [trm_ghost ghost_call]))

let trm_ghost_end (ghost_pair_var: var): trm =
  trm_add_cstyle GhostCall (trm_apps (trm_var ghost_end) [trm_var ghost_pair_var])

let trm_ghost_pair ?(name: string option) (ghost_call: ghost_call) : var * trm * trm =
  let ghost_pair_var = generate_ghost_pair_var ?name () in
  (ghost_pair_var,
   trm_ghost_begin ghost_pair_var ghost_call,
   trm_ghost_end ghost_pair_var)

let trm_ghost_custom_pair ?(name: string option) (forward_fn: trm) (backward_fn: trm): var * trm * trm =
  let ghost_pair_var = generate_ghost_pair_var ?name () in
  let ghost_with_reverse = trm_apps (trm_var with_reverse) [forward_fn; backward_fn] in
  (ghost_pair_var,
   trm_ghost_begin ghost_pair_var { ghost_fn = ghost_with_reverse; ghost_args = [] },
   trm_ghost_end ghost_pair_var)

let trm_ghost_scope ?(pair_name: string option) (ghost_call: ghost_call) (t: trm): trm =
  let _, ghost_begin, ghost_end = trm_ghost_pair ?name:pair_name ghost_call in
  Nobrace.trm_seq_nomarks [ghost_begin; t; ghost_end]

let trm_ghost_begin_inv (t: trm): (var * ghost_call) option =
  Pattern.pattern_match t [
    Pattern.(trm_let __ !__ __ (trm_apps1 (trm_var (var_eq ghost_begin)) (trm_apps !__ nil !__))) (fun pair_var ghost_fn ghost_args ->
      Some (pair_var, { ghost_fn; ghost_args })
    );
    Pattern.__ None
  ]

let trm_ghost_end_inv (t: trm): var option =
  Pattern.pattern_match t [
    Pattern.(trm_apps1 (trm_var (var_eq ghost_end)) (trm_var !__)) (fun pair_var -> Some pair_var);
    Pattern.__ None
  ]

let ghost_rewrite = toplevel_var "rewrite"

let trm_ghost_rewrite (before: formula) (after: formula) (by: formula): trm =
  trm_ghost (ghost_call ghost_rewrite ["H1", before; "H2", after; "by", by])

let ghost_forget_init = name_to_var "forget_init"

let trm_ghost_forget_init (f: formula): trm =
  trm_ghost (ghost_call ghost_forget_init ["H", f])
