
open Ast
open PPrint
open Target

(* debug flags *)
let debug = false
let debug_rec = false

(* mark used for marking trms that should be skipped by the simplifier *)
let mark_nosimpl = "__arith_core_nosimpl"

(* [has_mark_nosimplf t]: check if [t] should be skipped by the simplifier or not.*)
let has_mark_nosimpl (t : trm) : bool =
  trm_has_mark mark_nosimpl t
  (* LATER Ast.trm_has_mark *)

(* arithmetic operation type *)
type arith_op =
  | Arith_shift
  | Arith_scale

(* [transform_aux aop inv pre_cast post_cast u t]: shifts or scale the right hand
    side of a set operation with term [u]
    [aop] - a flag to decide if the arithmetic operation should be Arith_scale
       or Arith_shift
    [inv] - a flag for the sign(plus or minus) of shifting
    [u] - shift size
    [pre_cast] - casting of type [pre_cast] performed on the right hand side of the
      set operation before shifting
    [post_cast] - casting of type [post_cast] performed after shifting
    [t] - the ast of the set operation *)
let transform_aux (aop : arith_op) (inv : bool) (u : trm) (pre_cast : typ option)
  (post_cast : typ option)(t : trm) : trm =
  let binop_op = match aop with
    | Arith_shift -> if inv then Binop_sub else Binop_add
    | Arith_scale -> if inv then Binop_div else Binop_mul
    in
  let trm_apps_binop t1 t2 = trm_apps (trm_binop binop_op) [t1; t2] in
  match t.desc with
  | Trm_apps(f, [lhs; rhs]) when is_set_operation t ->
    begin match pre_cast, post_cast with
     | None, None -> trm_replace (Trm_apps (f, [lhs; trm_apps_binop rhs u])) t

     | None, Some ty -> trm_replace (Trm_apps (f, [lhs;trm_cast ty (trm_apps_binop rhs u)])) t

     | Some ty, None -> trm_replace (Trm_apps (f, [lhs;
                    trm_apps_binop (trm_cast ty rhs) u])) t
     | _ -> fail t.loc "Arith_core.transform_aux: can't apply both pre-casting
                        and post-casting"
    end
  | Trm_apps (_, [arg]) when is_get_operation t ->
    begin match pre_cast, post_cast with
     | None , None -> trm_apps_binop t u
     | None, Some ty -> trm_cast ty (trm_apps_binop t u)
     | Some ty, None -> trm_apps_binop (trm_cast ty t)  u
     | _ -> fail t.loc "Arith_core.transfom_aux: can't apply both pre-casting
                        and post-casting"
    end
  | _ -> fail t.loc "Arith_core.transform_aux: expected a get or a set operation"

(* [transform aop inv pre_cast post_cast u t p: applies [transform_aux] at the trm with path [p]. *)
let transform (aop : arith_op)(inv : bool) (u : trm) (pre_cast : typ option)
  (post_cast : typ option) : Transfo.local =
  apply_on_path (transform_aux aop inv u pre_cast post_cast)


(* [apply_aux op arg t]: applies binary_operation [op] on [t] with the second  argument of the operation being [arg],
    [op] -the binary operation to apply.
    [arg] - the second operand of [op].
    [t] - the first argument in the performed operation. *)
let apply_aux (op : binary_op) (arg : trm) (t : trm) : trm =
  trm_apps (trm_binop op) [t; arg]

(* [apply op arg t p]: applies [transform_aux] at the trm [t] with path [p]. *)
let apply (op : binary_op) (arg : trm) : Transfo.local =
  apply_on_path (apply_aux op arg)

(******************************************************************************)
(*                          Simplifier                                        *)
(******************************************************************************)

(* id type *)
type id = int

(* generate a new id *)
let next_id = Tools.fresh_generator ()

(* [expr]: expression type, it may be a literal expression, an atom expression or an arithmetic expression *)
type expr =
  | Expr_int of int
  | Expr_double of float
  | Expr_atom of id
  | Expr_sum of wexprs
  | Expr_prod of wexprs

(* [exprs]: list of expressions *)
and exprs = expr list

(* [wexprs]: weighted list of expressions *)
and wexprs = wexpr list

(* [wexpr]: weighted expression *)
and wexpr = (int * expr)

(* [Atom_map]: a map from atom ids to the corresponding terms *)
module Atom_map = Map.Make(Int)

(* [atom_map]: atom map for storing atoms *)
type atom_map = trm Atom_map.t

(* [no_atoms]: empty atom map *)
let no_atoms = Atom_map.empty

(******************************************************************************)
(*                          Smart constructors                                *)
(******************************************************************************)

(* [expr_mul we1 e]: multiplication expression  *)
let expr_mul (we1 : wexprs) (e : expr) : expr =
  match e with
  | Expr_prod [_,Expr_int 1] -> Expr_prod we1
  | Expr_prod wes -> Expr_prod (we1 @ wes)
  | _ -> Expr_prod  ((1,e) :: we1)

(* [expr_add we1 e]: addition expression  *)
let expr_add (we1 : wexprs) (e : expr) : expr =
  match e with
  | Expr_sum [_,Expr_int 0] -> Expr_sum we1
  | Expr_sum wes -> Expr_sum (we1 @ wes)
  | _ -> Expr_sum ((1, e) :: we1)

(* [expr_sum_nonweighted es]: nonweighted sum expression *)
let expr_sum_nonweighted (es : exprs) : expr =
   Expr_sum (List.map (fun e -> (1,e)) es)

(* [expr_prod_nonweighted es]: nonweighted prod expression *)
let expr_prod_nonweighted (es : exprs) : expr =
   Expr_prod (List.map (fun e -> (1,e)) es)


(* [apply_bottom_up]: is a combinator that takes a transformation and applies it recursively,
   bottom up through a term. *)
let rec apply_bottom_up (f : expr -> expr) (e : expr) : expr =
  let apply_wexprs (wes : wexprs) : wexprs =
    List.map (fun (w,e) -> (w, apply_bottom_up f e)) wes in
  match e with
  | Expr_sum wes -> f (Expr_sum (apply_wexprs wes))
  | Expr_prod wes -> f (Expr_prod (apply_wexprs wes))
  | _ -> f e


let identity (e : expr) : expr =
  e

(* [normalize_one e]:
   - collapses nested sums onto a single sum, and likewise for nested products
   - turns a product of an expression with a constant integer as a weighted
     expression in the parent sum
   - eliminates products and sums with a single expression of weight one
   - eliminates products and sums with an empty list
   - eliminates elements with weight zero
   - eliminates +0 is sums and *1 in produts *)
let normalize_one (e : expr) : expr =
  let e =
    match e with
    | Expr_sum wes ->
        Expr_sum (List.concat_map (function
          | (_ai, Expr_int 0) | (_ai, Expr_double 0.) -> []
          | (1, Expr_int n) -> [(n, Expr_int 1)]
          | (0, _ei) -> []
          | (1, Expr_sum wesi) -> wesi
          | (-1, Expr_sum [(ai,ei)]) -> [(-ai,ei)]
          | (ai, Expr_prod [(1, Expr_int bi); (1,ei)]) (* Optional? *)
          | (ai, Expr_prod [(1,ei); (1, Expr_int bi)]) (* Optional? *)
          | (ai, Expr_prod [(bi, Expr_int 1); (1,ei)])
          | (ai, Expr_prod [(1,ei); (bi, Expr_int 1)])  -> [(ai * bi, ei)]
          | we -> [we]) wes)
    | Expr_prod wes ->
        Expr_prod (List.concat_map (function
          | (_ai, Expr_int 1) | (_ai, Expr_double 1.) -> []
          | (0, _ei) -> []
          | (1, Expr_prod wesi) -> wesi
          | (-1, Expr_prod wesi) -> List.map (fun (w,ei) -> (-w, ei)) wesi
          | (ai, Expr_sum [(bi, Expr_int 1)]) -> [(ai, Expr_int bi)]
          | (ai, Expr_sum [(bi, ei)]) -> [(ai, Expr_int bi); (ai, ei)]
          | we -> [we]) wes)
    | _ -> e
    in
  match e with
  | Expr_sum [] -> Expr_int 0
  | Expr_prod [] -> Expr_int 1
  | Expr_sum [(1,e1)] -> e1
  | Expr_prod [(1,e1)] -> e1
  | _ -> e

(* [normaliez e]: applies [normalize_one] in a bottom up fashion *)
let normalize (e : expr) : expr =
  apply_bottom_up normalize_one e


(* [create_or_reuse_atom_for_trm atoms t]: auxiliary function for [trm_to_naive_expr]*)
let create_or_reuse_atom_for_trm (atoms : atom_map ref) (t : trm) : id =
  let no_id = -1 in
  let occ = ref no_id in
  Atom_map.iter (fun id tid ->
    if !occ = no_id && Internal.same_trm t tid then occ := id) !atoms;
    if !occ = no_id
      then begin
        let new_id = next_id() in
        atoms := Atom_map.add new_id t !atoms;
        occ := new_id
      end;
  !occ

(* [trm_to_naive_expr]: conversion of a trm from the AST into an expr, plus a map that for each atom gives
    the corresponding term *)
let trm_to_naive_expr (t : trm) : expr * atom_map =
  let atoms = ref Atom_map.empty in
  let rec aux (t : trm) : expr =
    let not_expression() = Expr_atom (create_or_reuse_atom_for_trm atoms t) in
    if has_mark_nosimpl t then not_expression() else
    match t.desc with
     | Trm_val (Val_lit (Lit_int n)) -> Expr_int n
     | Trm_val (Val_lit (Lit_double n)) -> Expr_double n
     | Trm_apps (f, [t1; t2]) ->
       begin match trm_prim_inv f with
        | Some (Prim_binop b) ->
          begin match b with
          | Binop_add | Binop_sub ->
             let w = match b with | Binop_add -> 1 | Binop_sub -> -1 | _ -> assert false in
             Expr_sum [(1, aux t1); (w, aux t2)]
          | Binop_mul | Binop_div ->
             let w = match b with | Binop_mul -> 1 | Binop_div -> -1 | _ -> assert false in
             Expr_prod [(1, aux t1); (w, aux t2)]
          | _ -> not_expression()
          end
        | _ -> not_expression()
       end
     | Trm_apps (f, [t1]) ->
       begin match trm_prim_inv f with
        | Some (Prim_unop Unop_neg) -> Expr_sum [(-1, aux t1)]
        | _ -> not_expression()
       end
     | _ -> not_expression()
     in
    let res = aux t in
    res, !atoms

(* [is_one e]: checks if e == 1 *)
let is_one (e : expr) : bool =
  match e with
  | Expr_int 1 | Expr_double 1.0 -> true
  | _ -> false

(* [parens_if_neg n d]: if [n] is negative then it add parentheses around [d] *)
let parens_if_neg (n:int) (d:document) : document =
  if n < 0 then parens d else d

(* [expr_to_string atoms e]: convert an expression to a string, in AST form *)
let expr_to_string (atoms : atom_map) (e : expr) : string =
  let rec aux (e : expr) : document =
    let auxw ((w,e) : wexpr) : document =
      parens ((string (string_of_int w)) ^^ comma ^^ aux e) in
    let auxwes (we : wexprs) : document =
      Tools.list_to_doc (List.map auxw we) in
    match e with
    | Expr_int n -> string (string_of_int n)
    | Expr_double n -> string (string_of_float n)
    | Expr_sum wes -> string "Sum" ^^ (auxwes wes)
    | Expr_prod wes -> string "Prod" ^^ (auxwes wes)
    | Expr_atom id ->
        begin match Atom_map.find_opt id atoms with
        | Some t1 ->
            begin match t1.desc with
            | Trm_var (_, x) -> string x.qvar_var
            | Trm_apps ({desc = Trm_val (Val_prim (Prim_unop Unop_get)); _},
               [{desc = Trm_var (_, x); _}]) -> string x.qvar_var
            | _ -> braces (AstC_to_c.trm_to_doc t1)
            end
        | _  ->
          if id > 26*26 then braces (string (string_of_int id)) else
            let string_of_small_int i =
              let c = Char.chr ((Char.code 'a') + i) in
              String.make 1 c in
            let s1 = (if id >= 26 then string_of_small_int (id / 26) else "") in
            let s2 = string_of_small_int (id mod 26) in
            braces (string "!" ^^ string (s1 ^ s2))
        end
  in
  Tools.document_to_string (aux e)


(* [expr_to_math_string atoms e]: converts an expression to a string, using mathematical notations *)
let expr_to_math_string (atoms : atom_map) (e : expr) : string =
  let power_to_doc (base : document) (power : int) : document =
     base ^^ string "^" ^^ string (string_of_int power)
     in
  let rec aux (e : expr) : document =
    match e with
    | Expr_int n -> string (string_of_int n)
    | Expr_double n -> string (string_of_float n)
    | Expr_sum we ->
      begin match we with
      | [] -> Printf.printf "Arith_core.expr: Expr_sum [] should never appear";
        (string (string_of_int 0))
      | _ ->
        let we_l = List.map (fun (w, e) ->
          if is_one e then
            string (string_of_int w)
          else begin
            let s = aux e in
            if w = 1
              then s
              else ((parens_if_neg w (string (string_of_int w))) ^^ star ^^ s)
          end
        ) we in
        Tools.list_to_doc ~sep:plus ~bounds:[lparen; rparen] we_l
      end
    | Expr_prod we ->
      begin match we with
      | [] -> Printf.printf "Arith_core.expr: Expr_prod [] should never appear";
        string (string_of_int 1)
      | _ ->
        let we_l = List.map (fun (w, e) ->
        power_to_doc (aux e) w
      ) we in
      Tools.list_to_doc ~sep:star ~bounds:[lparen; rparen] we_l
      end
    | Expr_atom id ->
      begin match Atom_map.find_opt id atoms with
      | Some t1 -> (AstC_to_c.trm_to_doc t1)
      | _  -> fail None "Arith_core.expr_to_math_string: couldn't convert
                        an atom expr to a trm"
      end
  in
  Tools.document_to_string (aux e)


(* [trm_to_expr t]: convert trm [t] to an expression*)
let trm_to_expr (t : trm) : expr * atom_map =
  let expr, atoms = trm_to_naive_expr t in
  if debug
    then Printf.printf "Expr after conversion: %s\n" (expr_to_string atoms expr);
  let res = normalize expr in
  if debug
    then Printf.printf "Expr after normalization: %s\n" (expr_to_string atoms res);
  res, atoms

(* [expr_to_trm atoms e]: converts expr [e] to trm  *)
let expr_to_trm (atoms : atom_map) (e : expr) : trm =
  let rec aux (e : expr) : trm =
    match e with
    | Expr_int n -> trm_int n
    | Expr_double n -> trm_double n
    | Expr_sum we ->
      if we = [] then failwith "expr_to_trm: assumes a normalized term";
      Xlist.fold_lefti (fun i acc (w, e) ->
        let (w, e) = match (w,e) with (n,Expr_int 1) -> (1,Expr_int n) | _ -> (w,e) in
        let c = if i = 0 then w else abs w in
        let y = aux e in
        let x =
          if c = 1
            then y
            else if c = -1 then trm_apps (trm_unop Unop_minus) [y]
            else trm_apps (trm_binop Binop_mul) [trm_int c; y]
          in
        if i = 0
          then x
          else trm_apps (trm_binop (if w >= 0 then Binop_add else Binop_sub)) [acc; x]
      ) (trm_unit ()) we
    | Expr_prod  we ->
      if we = [] then failwith "expr_to_trm: assumes a normalized term";
      let rec power t n =
        if n = 0
          then failwith "Arith_core.expr_to_trm: assumes a normalized term"
          else if n < 0 then trm_apps (trm_binop Binop_div) [trm_double 1.0; power t (-n)]
          else if n = 1 then t
          else trm_apps (trm_binop Binop_mul) [t; power t (n-1)]
        in
      Xlist.fold_lefti (fun i acc (w,e) ->
        if i = 0
          then power (aux e) w
          else trm_apps (trm_binop (if w > 0
             then Binop_mul
             else Binop_div)) [acc; power (aux e) (abs w)]
      ) (trm_unit ()) we
    | Expr_atom id ->
        begin match Atom_map.find_opt id atoms with
        | Some t1 -> t1
        | _ -> fail None "Arith_core.expr_to_trm: couldn't convert an atom expr to a trm"
        end
    in
  aux e


(* [cleanup_true]: perform cleanup *)
let cleanup_true = true

(* [cleanup_false]: don't perform cleanup *)
let cleanup_false = false


(* [apply_bottom_up_if]: is a combinator for either applying a transformation recursively
   or applying it only at the top level, according to the [recurse] argument.
   If the [cleanup] argument is true, then after each call to the transformation,
   the operation [normalize_one] is called. *)
let apply_bottom_up_if (recurse : bool) (cleanup : bool) (f : expr -> expr)
  (e : expr) : expr =
  let f_with_cleanup e =
    let e1 = (if cleanup then normalize_one else identity) e in
    let e2 = f e1 in
    let e3 = (if cleanup then normalize_one else identity) e2 in
    if debug_rec
      then Printf.printf "Step:\n\t%s\n\t%s\n\t%s\n" (expr_to_string no_atoms e)
        (expr_to_string no_atoms e2) (expr_to_string no_atoms e3);
    e3 in
  if recurse
    then apply_bottom_up f_with_cleanup e
    else f_with_cleanup e

(* [apply_bottom_up_debug e]: function used only for debugging purposes *)
let apply_bottom_up_debug (e : expr) : expr =
  let f ei =
    if debug then Printf.printf "Bottom-up %s\n" (expr_to_string no_atoms ei);
    ei in
  apply_bottom_up f e

(* [apply_bottom_up_if_debug]: function used only for debugging purposes *)
let apply_bottom_up_if_debug (recurse : bool) (cleanup : bool) (e : expr) : expr =
  let f ei =
    let ej = (if cleanup then normalize_one else identity) ei in
    if debug
      then Printf.printf "Bottom-up-if:\n\t%s\n\t%s\n"
            (expr_to_string no_atoms ei) (expr_to_string no_atoms ej); ej
    in
  apply_bottom_up_if recurse cleanup f e

(* LATER: Use a map instead of a list *)
(* [gather_one e]: regroups similar expression that appear inside a same product
    or sum. For example, [2 * e1 + (-1)*e1] simplifies to [e1] and
    [e1 * e2 * e1^(-1)] simplifies to [e2]. *)
let gather_one (e : expr) : expr =
  let rec insert (acc : wexprs) ((w,e) : wexpr) : wexprs =
      match acc with
      | [] -> [(w,e)]
      | (wi,ei)::acc2 -> if e = ei
          then (wi+w,ei)::acc2
          else (wi,ei)::(insert acc2 (w,e))
    in
  let gather_wexprs (wes : wexprs) : wexprs =
     List.fold_left insert [] wes
    in
  match e with
  | Expr_sum wes -> Expr_sum (gather_wexprs wes)
  | Expr_prod wes -> Expr_prod (gather_wexprs wes)
  | _ -> e

(* [gather_common recurse_bool e]: apply gather one in a full expression
    if recurse is set to true *)
let gather_common (recurse : bool) (e : expr) : expr =
  apply_bottom_up_if recurse cleanup_true gather_one e

let gather = gather_common false
let gather_rec = gather_common true

(* [expand_one e]: expands sums that appear inside product.
    For example, [e1 * (e2 + e3)] becomes [e1 * e2 + e1 * e3]
    The function is equeal to the identity if no expansion can be performed.
    At the very end it will apply [normalize] to the result.*)
let expand_one (e : expr) : expr =
  (* [acc] corresponds to the list of terms in the current sum;
     [e^w] is the term to distribute over the sum described by [acc].
     For example, assume [acc] is [e1; e2].
     If [e^w = w1*a1 + w2*a2], then we produce
      [w1*a1 * e1 + w1*a1 * e2 + w2*a2 * e1 + w2*a2 * e2].
     Else, we produce [e^w*e1; e^w*e2]. *)
  let aux ((w,e) : wexpr) (acc : exprs) : exprs =
    match (w,e) with
    | 1, (Expr_sum wes) ->
        List.concat_map (fun (wk,ak) ->
          List.map (fun ei ->
            expr_prod_nonweighted [(Expr_int wk); ak; ei]) acc
        ) wes
    | _ ->
      List.map (fun ei -> expr_mul [(w,e)] ei) acc
    in
  let r = match e with
    | Expr_prod wes ->
        let exprs_in_sum = List.fold_right aux wes [Expr_int 1] in
        expr_sum_nonweighted exprs_in_sum
    | _ -> e
    in
  normalize r

(* [expand_common recurse e]: calls [expand_one] recursively, calling the [gather] operations
    after each step. *)
let expand_common (recurse : bool) (e : expr) : expr =
  let tr (ei : expr) : expr =
    gather_rec (expand_one ei) in
  apply_bottom_up_if recurse cleanup_true tr e

let expand = expand_common false
let expand_rec = expand_common true (* Warning: quadratic, because normalize
                                        all and gather_rec at each step *)

(* [map_on_arith_nodes tr t]: applies arithmetic simplification [tr] in depth of [t]*)
let rec map_on_arith_nodes (tr : trm -> trm) (t : trm) : trm =
  if has_mark_nosimpl t
    then t
    else if is_prim_arith_call t then tr t
    else
      trm_map (map_on_arith_nodes tr) t


(* DEBUG let c = ref 0 *)
let simplify_at_node (f_atom : trm -> trm) (f : expr -> expr) (t : trm) : trm =
  let expr, atoms = trm_to_expr t in

  let atoms2 =
    (* If we could not extract any structure, then we don't process recursively the atoms
       using [f_atom], else we would trigger an infinite loop when [indepth=true]. *)
    match expr with
    | Expr_atom _id -> atoms
    | _ ->  Atom_map.map f_atom atoms
    in
  let expr2 = f expr in
  if debug then Printf.printf "Expr after transformation: %s\n" (expr_to_string atoms expr2);
  expr_to_trm atoms2 expr2

(* [simplify_aux indepth f t]: converts node [t] to an expression, then applies the
     simplifier [f], then it converts it back to a trm
    params:
      [f]: simplifier function
      [t]: the node on which the simplifications should be performed
    return:
      update t with the simplified expressions
  LATER: should [simplify_aux false f t] fail if [t] is not an application of prim_arith? *)
let rec simplify_aux (indepth : bool) (f : expr -> expr) (t : trm) : trm =
  if not indepth
    then begin
       let f_atom_identity = (fun ti -> ti) in
       simplify_at_node f_atom_identity f t
      end
    else begin
     let f_atom_simplify = simplify_aux indepth f in
     map_on_arith_nodes (simplify_at_node f_atom_simplify f) t end

(* [simplify indepth f t p]: applies [simplify_aux] at the trm with path [p] *)
let simplify (indepth : bool) (f : expr -> expr) : Transfo.local =
  apply_on_path (simplify_aux indepth f)
