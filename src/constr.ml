open Ast
open Str
open Tools
open Path


(******************************************************************************)
(*                        Data structure for targets                          *)
(******************************************************************************)

(* [trm_kind]: type to classify trms into four main classes
    1) Structuring statements
    2) Instructions
    3) Expression
    4) others  *)
type trm_kind =
  | TrmKind_Typedef (* type definition that appears in the AST *)
  | TrmKind_Ctrl    (* control statement, for loops, while loops and unit-type if statements  *)
  | TrmKind_Instr
  | TrmKind_Expr
  | TrmKind_Any



(* [rexp]: data structure for storing regular expressions *)
type rexp = {
  rexp_desc : string; (* printable version of regexp *)
  rexp_exp : regexp;
  rexp_substr : bool;
  rexp_trm_kind : trm_kind }

(* [target_relative]: target kind *)
type target_relative =
    | TargetAt
    | TargetFirst
    | TargetLast
    | TargetBefore
    | TargetAfter

(* [target_occurrences]: the number of targets to match *)
type target_occurrences =
    | ExpectedOne  (* = 1 occurence, the default value *)
    | ExpectedNb of int  (* exactly n occurrences *)
    | ExpectedMulti  (* > 0 number of occurrences *)
    | ExpectedAnyNb  (* any number of occurrences *)
    | ExpectedSelected of int option * int list
    | FirstOcc
    | LastOcc
(* A [target] is a list of constraints to identify nodes of the AST
   that we require the result path to go through. *)

(* [target]: is a list of constraints to identify nodes of the AST that we require the result path to go through. *)
type target = constr list

(* [targets]: a list of targets *)
and targets = target list

(* [depth]: targets by default resolve at any depth, but they can also be resolved only on a specific depth *)
and depth =
  | DepthAny
  | DepthAt of int

(* [constr]: are the unit constrainst that are checked when target is being resolved, a constraint can be one of the following ones
     - direction constraints
     - include directives constraints
     - depth constraints
     - regexp constraints
     - trm constraints
     - type constraints
     - logic constraints
     - occurrence constraints
     - relative constraints
     - OpenMP pragam constraints *)
and constr =
  | Constr_depth of depth
  | Constr_dir of dir
  | Constr_include of string
  | Constr_regexp of rexp
  (* for: init, cond, step, body *)
  | Constr_for_c of target * target * target * target
  (* for index, start, stop, step, body *)
  | Constr_for of constr_name * target * loop_dir option * target * target * target
  (* while: cond, body *)
  | Constr_while of target * target
  (* do while: body, cond *)
  | Constr_do_while of target * target
  (* if: cond, then, else *)
  | Constr_if of target * target * target
  (* decl_var: name, body *)
  | Constr_decl_var of typ_constraint * constr_name * target
  (* decl_fun: name, args, body is_def*)
  | Constr_decl_fun of typ_constraint * constr_name * target_list_pred * target * bool
  (* decl_type: name *)
  | Constr_decl_type of constr_name
  (* decl_enum: name, constants *)
  | Constr_decl_enum of constr_name * constr_enum_const
  (* seq *)
  | Constr_seq of target_list_pred
  (* var *)
  | Constr_var of constr_name
  (* lit *)
  | Constr_lit of (lit -> bool)
  (* app: function, arguments *)
  | Constr_app of target * target_list_pred * bool
  (* label *)
  | Constr_label of constr_name * target
  (* goto *)
  | Constr_goto of constr_name
  (* return *)
  | Constr_return of target
  (* abort *)
  | Constr_abort of abort_kind
  (* accesses: base, accesses *)
  | Constr_access of target * constr_accesses * bool
  (* switch: cond, cases *)
  | Constr_switch of target * constr_cases
  (* Target relative to another trm *)
  | Constr_relative of target_relative
  (* Number of  occurrences expected  *)
  | Constr_occurrences of target_occurrences
  (* List of constraints *)
  | Constr_target of constr list
  (* Constraint used for argument match *)
  | Constr_bool of bool
  (* Constraint that matches only the root of the AST *)
  | Constr_root
  (* Constraint that matches primitive operations *)
  | Constr_prim of (prim -> bool)
  (* Constraint that matches ast nodes whose marks satisfy the predicate *)
  | Constr_mark of (mark -> bool) * string
  (* Constraint that matches a union of targets *)
  | Constr_or of target list
  (* Constraint that matches an intersection of targets *)
  | Constr_and of target list
  (* Constrain that match the diff between two targets *)
  | Constr_diff of target list * target list
  (* Constraint to match arguments  that sastify both the name predicated and the type predicate *)
  | Constr_arg of var_constraint * typ_constraint
  (* Constraint to match ast nodes of types that satisfy the type predicate *)
  | Constr_hastype of typ_constraint
  (* Constraint to match variable initialization value *)
  | Constr_var_init
  (* Constraint to match an array initialization list *)
  | Constr_array_init
  (* Constraint to match an struct initialization list  *)
  | Constr_struct_init
  (* Constraint to match an omp directive *)
  | Constr_omp of (directive->bool) * string

(* LATER: optimize constr_of_path; should be recognized by resolution,
   and processed more efficiently; checking that the start of the path
   is the root, then reaching directly the position, assuming the path
   is valid (option: raise an exception if the path is invalid). *)

(* [typ_constraint]: predicate on types *)
and typ_constraint = typ -> bool

(* [arg_constraint]: constraint on an argument, represented as a target with a single item
    of the form [cArg ..] or [cTrue]  *)
and arg_constraint = target

(* [var_constraint]: constraint over variable names *)
and var_constraint = string -> bool

(* [constr_name): names involved in constraints, e.g. for goto labels *)
and constr_name = rexp option

(* [constr_enum_const]: constraint on const enums *)
and constr_enum_const = ((constr_name * target) list) option

(* [constr_accesses]: constraint of a list of struct accesses or array accesses *)
and constr_accesses = (constr_access list) option

(* Predicate for expressing constraints over a list of subterms. It consists of:
   - a function that produces the constraints for the i-th subterm
   - a function that takes a list of booleans indicating which of the subterms
     matched their respective constraint, and returns a boolean indicating whether
     the full list should be considered as matching or not.
   - a string that explains what was the user intention *)

(* [target_list_pred]: predicate for expressing constraints over a list of stubterms. It consists of
    - a function that produces the constraints for the i-th subterm
    - a function that takes a list of booleans indicating which of the subterms match their respective constraint,
        and returns a boolean indicatin whether the full list should be considered as matching or not
    - a string that explains what was the user intention *)
and target_list_pred =
  { target_list_pred_ith_target : int -> target;
    target_list_pred_validate : bool list -> bool;
    target_list_pred_to_string : unit -> string; }

(* [constr_access]: constraint over a single struct or array access  *)
and constr_access =
  (* array indices may be arbitrary terms *)
  | Array_access of target
  (* struct fields are strings *)
  | Struct_access of constr_name
  | Any_access

(* [constr_cases]: for each case, its kind and a constraint on its body *)
and constr_cases = ((case_kind * target) list) option

(* [case_kind]: switch case kind used on [constr_cases] *)
and case_kind =
  (* case: value *)
  | Case_val of target
  | Case_default
  | Case_any

(* [abort_kind]: abort instruction kind used for abort insitruction constraints  *)
and abort_kind =
  | Any
  | Return
  | Break
  | Continue


(* [target_simple]: is a [target] without any of the following relative constraints
   Constr_relative, Constr_occurrences, Constr_target.
   Note: It can include [cStrict] *)

type target_simple = target



(* [target_struct]: structured representation of a [target] that decomposes the special constructors
   such as Constr_relative, Constr_occurrences, Constr_target from the [target_simple]. *)
type target_struct = {
   target_path : target_simple; (* this path contains no nbMulti/nbEx/tBefore/etc.., only cStrict can be there *)
   target_relative : target_relative;
   target_occurrences : target_occurrences; }


(* [make_target_list_pred ith_target validate to_string]: creates a simple target_pred object with with its components
   being [ith_target], [validate] and [to_string] *)
let make_target_list_pred (ith_target : int -> target) (validate : bool list -> bool) (to_string : unit -> string) : target_list_pred =
  { target_list_pred_ith_target = ith_target;
    target_list_pred_validate = validate;
    target_list_pred_to_string = to_string; }

(* [depth_pred d]: predicate on the depth *)
let depth_pred (d : depth) : depth =
  match d with
  | DepthAny -> DepthAny
  | DepthAt n ->
      if n <= 0 then fail None "Constr.depth_pred: argument of DepthAt is not positive";
      DepthAt (n-1)

(******************************************************************************)
(*                        Pretty-printing of targets                          *)
(******************************************************************************)

(* [trm_kind_to_string k]: pretty prints trm kind [k] *)
let trm_kind_to_string (k : trm_kind) : string =
  match k with
  | TrmKind_Typedef -> "Typedef"
  | TrmKind_Ctrl -> "Ctrl"
  | TrmKind_Instr -> "Instr"
  | TrmKind_Expr -> "Expr"
  | TrmKind_Any -> "Any"

(* [rexp_to_string r]: pretty prints rexp [r] *)
let rexp_to_string (r : rexp) : string =
  (trm_kind_to_string r.rexp_trm_kind) ^ "-" ^
  (if r.rexp_substr then "Sub" else "Exact") ^ "-" ^
  r.rexp_desc

(* [depth_to_string depth]: pretty prints depth [depth] *)
let depth_to_string (depth : depth) : string =
  match depth with
  | DepthAny -> "DepthAny"
  | DepthAt n -> "DepthAt " ^ string_of_int n

(* [constr_to_string c]: pretty prints constraint [c] *)
let rec constr_to_string (c : constr) : string =
  match c with
  | Constr_depth depth -> "Depth " ^ (depth_to_string depth)
  | Constr_dir d -> dir_to_string d
  | Constr_include s -> "Include " ^ s
  | Constr_regexp r -> "Regexp " ^ rexp_to_string r
  | Constr_for_c (p_init, p_cond, p_step, p_body) ->
     let s_init = target_to_string p_init in
     let s_cond = target_to_string p_cond in
     let s_step = target_to_string p_step in
     let s_body = target_to_string p_body in
     "For (" ^ s_init ^ ", " ^ s_cond ^ ", " ^ s_step ^ ", " ^ s_body ^ ")"
  | Constr_for (p_index, p_start, p_direction, p_stop, p_step, p_body) ->
    let s_index =
      match p_index with | None -> "_" | Some r -> rexp_to_string r
    in
    let s_direction = match p_direction with
    | Some DirUp -> "Up"
    | Some DirUpEq -> "UpEq"
    | Some DirDown -> "Down"
    | Some DirDownEq -> "DownEq"
    | None -> "AnyDirection"
    in
    let s_start = target_to_string p_start in
    let s_stop = target_to_string p_stop in
    let s_step = target_to_string p_step in
    let s_body = target_to_string p_body in
    "For ("^s_index ^ " ," ^ s_direction ^ " , " ^ s_start ^ ", " ^ s_stop ^ ", " ^ s_step ^ ", " ^ s_body ^ ")"
  | Constr_while (p_cond, p_body) ->
     let s_cond = target_to_string p_cond in
     let s_body = target_to_string p_body in
     "While (" ^ s_cond ^ ", " ^ s_body ^ ")"
  | Constr_do_while (p_body, p_cond) ->
     let s_body = target_to_string p_body in
     let s_cond = target_to_string p_cond in
     "Do_While (" ^ s_body ^ ", " ^ s_cond ^ ")"
  | Constr_if (p_cond, p_then, p_else) ->
     let s_cond = target_to_string p_cond in
     let s_then = target_to_string p_then in
     let s_else = target_to_string p_else in
     "If (" ^ s_cond ^ ", " ^ s_then ^ ", " ^ s_else ^ ")"
  | Constr_decl_var (_ty_pred, name, p_body) ->
     let s_name =
       match name with | None -> "_" | Some r -> rexp_to_string r
     in
     let s_body = target_to_string p_body in
     "Decl_var (<ty_pred>, " ^ s_name ^ ", " ^ s_body ^ ")"
     (* LATER: add a string representation for type constraints *)
  | Constr_decl_fun (_ty_pred,name, _tgt_list_pred, p_body, is_def) ->
    let s_name =
       match name with | None -> "_" | Some r -> rexp_to_string r
     in
     let spred = _tgt_list_pred.target_list_pred_to_string() in
     let s_body = target_to_string p_body in
     let s_is_def = string_of_bool is_def in
     "Decl_fun (<ty_pred>, " ^ s_name ^ spred ^ s_body ^ s_is_def ^ ")"

  | Constr_decl_type name ->
     let s_name =
       match name with | None -> "_" | Some r -> rexp_to_string r
     in
     "Decl_type " ^ s_name
  | Constr_decl_enum (name, c_const) ->
     let s_name =
       match name with | None -> "_" | Some r -> rexp_to_string r
     in
     let s_const =
       match c_const with
       | None -> "_"
       | Some np_l ->
          let sl =
            List.map
              (fun (n, p) ->
                let s_n =
                  match n with | None -> "_" | Some r -> rexp_to_string r
                in
                let s_p = target_to_string p in
                "(" ^ s_n ^ ", " ^ s_p ^ ")"
              )
              np_l
          in
          list_to_string sl
     in
     "Decl_enum (" ^ s_name ^ ", " ^ s_const ^ ")"
  | Constr_seq tgt_list_pred ->
     let spred = tgt_list_pred.target_list_pred_to_string() in
     Printf.sprintf "Seq (%s)" spred
  | Constr_var name ->
     "Var " ^ (match name with | None -> "_" | Some r -> rexp_to_string r)
  | Constr_lit _ -> "Lit"
  | Constr_app (p_fun, tgt_list_pred, accept_encoded) ->
    let spred = tgt_list_pred.target_list_pred_to_string() in
    let s_fun = target_to_string p_fun in
    "App (" ^ s_fun ^ ", " ^ spred ^ string_of_bool(accept_encoded) ^ ")"
  | Constr_label (so, p_body) ->
     let s_label =
       match so with | None -> "_" | Some r -> rexp_to_string r
     in
     let s_body = target_to_string p_body in
     "Label (" ^ s_label ^ ", " ^ s_body ^ ")"
  | Constr_goto so ->
     let s_label =
       match so with | None -> "_" | Some r -> rexp_to_string r
     in
     "Goto " ^ s_label
  | Constr_return p_res ->
      let s_res = target_to_string p_res in
      "Return " ^ s_res
  | Constr_abort kind ->
     let s_kind =
       match kind with
       | Any -> "Any"
       | Return -> "Return"
       | Break -> "Break"
       | Continue -> "Continue"
     in
     "Abort_" ^ s_kind
  | Constr_access (p_base, ca, _) ->
     let s_accesses =
       match ca with
       | None -> "_"
       | Some cal -> list_to_string (List.map access_to_string cal)
     in
     let s_base = target_to_string p_base in
     "Access (" ^ s_accesses ^ ", " ^ s_base ^ ")"
  | Constr_switch (p_cond, cc) ->
     let s_cond = target_to_string p_cond in
     let s_cases =
       match cc with
       | None -> "_"
       | Some cl ->
          let string_of_kind = function
            | Case_val p_val -> "Case_val " ^ target_to_string p_val
            | Case_default -> "Default"
            | Case_any -> "Any"
          in
          let sl =
            List.map
              (fun (k, p_body) ->
                let s_body = target_to_string p_body in
                "(" ^ string_of_kind k ^ ", " ^ s_body ^ ")"
              )
              cl
          in
          list_to_string sl
     in
     "Switch (" ^ s_cond ^ ", " ^ s_cases ^ ")"
  | Constr_relative tr -> target_relative_to_string tr
  | Constr_occurrences oc -> target_occurrences_to_string oc
  | Constr_target cl ->
    let string_cl = List.map constr_to_string cl in
    "Target" ^ list_to_string string_cl
  | Constr_bool b -> if b then "True" else "False"
  | Constr_root -> "Root"
  | Constr_prim _ -> "Prim"
  | Constr_mark (_, str) -> "Mark (" ^ str ^ ")"
  | Constr_or tl -> "Or (" ^ Tools.list_to_string (List.map target_to_string tl) ^ ")"
  | Constr_and tl -> "And (" ^ Tools.list_to_string (List.map target_to_string tl) ^ ")"
  | Constr_diff (tl1, tl2) -> "Diff (" ^ Tools.list_to_string (List.map target_to_string tl1) ^ "," ^ Tools.list_to_string (List.map target_to_string tl2) ^ ")"
  | Constr_arg _ -> "<Constr_args>"
  | Constr_hastype _ -> "<Constr_hastype>"
  | Constr_var_init -> "Var_init"
  | Constr_array_init -> "Array_init "
  | Constr_struct_init -> "Struct_init"
  | Constr_omp (_, str) -> "Omp (" ^ str ^ ")"



(* [target_to_string tg]: pretty prints target [tg] *)
and target_to_string (tg : target) : string =
  list_to_string (List.map constr_to_string tg)


(* [target_struct_to_string tgs]: pretty prints target struct [tgs] *)
and target_struct_to_string (tgs : target_struct) : string =
  "TargetStruct(" ^
    target_relative_to_string tgs.target_relative ^ ", " ^
    target_occurrences_to_string tgs.target_occurrences ^ ", " ^
    target_to_string tgs.target_path ^ ")"

(* [target_occurrences_to_string occ]: pretty prints target_occurrences [occ] *)
and target_occurrences_to_string (occ : target_occurrences) =
  match occ with
  | FirstOcc -> "FirstOcc"
  | LastOcc -> "LastOcc"
  | ExpectedOne -> "ExpectedOne"
  | ExpectedNb n -> Printf.sprintf "ExpectedNb(%d)" n
  | ExpectedMulti -> "ExpectedMulti"
  | ExpectedAnyNb -> "ExpectedAnyNb"
  | ExpectedSelected (ex_opt, il) ->
    let exact_nb_s = match ex_opt with
    | None -> "None"
    | Some i -> "Some " ^ (string_of_int i) in
    Printf.sprintf "ExpectedSelect(%s, %s)" exact_nb_s (Tools.list_to_string (List.map (string_of_int) il))

(* [target_relative_to_string rel]: pretty prints the relative target [rel] *)
and target_relative_to_string (rel : target_relative) =
  match rel with
  | TargetAt -> "TargetAt"
  | TargetFirst -> "TargetFirst"
  | TargetLast -> "TargetLast"
  | TargetBefore -> "TargetBefore"
  | TargetAfter -> "TargetAfter"

(* [access_to_string ca]: pretty prints access [ca] *)
and access_to_string (ca : constr_access) : string =
  match ca with
  | Array_access p_index ->
     let s_index = target_to_string p_index in
     "Array_access " ^ s_index
  | Struct_access so ->
     let s_field =
       match so with | None -> "_" | Some r -> rexp_to_string r
     in
     "Struct_access " ^ s_field
  | Any_access -> "Any_access"

(* [constr_map f c]: applies [f] recursively on constraint [c] *)
let constr_map (f : constr -> constr) (c : constr) : constr =
  (* LATER: optimize using identity reconstruction, like trm_map with optim *)
  let aux (cs:target) : target =
    List.map f cs in
  let auxs (tgs:targets) : targets =
    List.map aux tgs in
  match c with
  | Constr_depth _depth -> c
  | Constr_dir _d -> c
  | Constr_include _s -> c
  | Constr_regexp _r -> c
  | Constr_for_c (p_init, p_cond, p_step, p_body) ->
     let s_init = aux p_init in
     let s_cond = aux p_cond in
     let s_step = aux p_step in
     let s_body = aux p_body in
     Constr_for_c (s_init, s_cond, s_step, s_body)
  | Constr_for (p_index, p_start, p_direction, p_stop, p_step, p_body) ->
      let s_start = aux p_start in
      let s_stop = aux p_stop in
      let s_step = aux p_step in
      let s_body = aux p_body in
      Constr_for (p_index, s_start, p_direction, s_stop, s_step, s_body)
  | Constr_while (p_cond, p_body) ->
     let s_cond = aux p_cond in
     let s_body = aux p_body in
     Constr_while (s_cond, s_body)
  | Constr_do_while (p_body, p_cond) ->
     let s_body = aux p_body in
     let s_cond = aux p_cond in
     Constr_do_while (s_body, s_cond)
  | Constr_if (p_cond, p_then, p_else) ->
     let s_cond = aux p_cond in
     let s_then = aux p_then in
     let s_else = aux p_else in
     Constr_if (s_cond, s_then, s_else)
  | Constr_decl_var (ty_pred, name, p_body) ->
     let s_body = aux p_body in
     Constr_decl_var (ty_pred, name, s_body)
  | Constr_decl_fun (ty_pred,name, tgt_list_pred, p_body, is_def) ->
     let s_body = if is_def then aux p_body else p_body in
     Constr_decl_fun (ty_pred,name, tgt_list_pred, s_body, is_def)
  | Constr_decl_type name -> c
  | Constr_decl_enum (name, c_const) ->
     let s_const = Tools.option_map (List.map (fun (n,tg) -> (n,aux tg))) c_const in
     Constr_decl_enum (name, s_const)
  | Constr_seq _tgt_list_pred -> c
  | Constr_var _name -> c
  | Constr_lit _-> c
  | Constr_app (p_fun, tgt_list_pred, accept_encoded) ->
     let s_fun = aux p_fun in
     Constr_app (s_fun, tgt_list_pred, accept_encoded)
  | Constr_label (so, p_body) ->
     let s_body = aux p_body in
     Constr_label (so, s_body)
  | Constr_goto so -> c
  | Constr_return p_res ->
      let s_res = aux p_res in
      Constr_return s_res
  | Constr_abort kind -> c
  | Constr_access (p_base, ca, inner) ->
     let s_base = aux p_base in
     Constr_access (s_base, ca, inner)
  | Constr_switch (p_cond, cc) ->
     let s_cond = aux p_cond in
     let s_cc = Tools.option_map (List.map (fun (k,tg) -> (k,aux tg))) cc in
     Constr_switch (s_cond, s_cc)
  | Constr_relative tr -> c
  | Constr_occurrences oc -> c
  | Constr_target cl ->
      Constr_target (aux cl)
  | Constr_bool b -> c
  | Constr_root -> c
  | Constr_prim _ -> c
  | Constr_mark (_, str) -> c
  | Constr_or tl ->
      Constr_or (auxs tl)
  | Constr_and tl ->
      Constr_and (auxs tl)
  | Constr_diff (tl1, tl2) ->
      let s_tl1 = auxs tl1 in
      let s_tl2 = auxs tl2 in
      Constr_diff (s_tl1, s_tl2)
  | Constr_arg _ -> c
  | Constr_hastype _ -> c
  | Constr_var_init -> c
  | Constr_array_init -> c
  | Constr_struct_init -> c
  | Constr_omp _ -> c

(* [get_target_regexp_kinds tgs]: gets the list of trm_kinds of the terms
   for which we would potentially need to use the string representation,
   for resolving the targets [tgs]. The result is either a list of kinds
   without [TrmKind_Any], or a singleton list made of [TrmKind_Any]. *)
let get_target_regexp_kinds (tgs : target list) : trm_kind list =

  (* LATER: could be optimize by avoiding the construction of a copy of c;
     but this should be done automatically when constr_map is optimized *)
  (* Note: we use lists, but this is fine because the lists are very short *)
  let res = ref [] in
  let add (k : trm_kind) : unit =
    match k with
    | TrmKind_Any -> res := [TrmKind_Any]
    | _ -> if not (List.mem k !res) then res := k :: !res
    in
  let rec explore (c : constr) : constr = (* LATER: optimize using a constr_iter instead of constr_map *)
    begin match c with
    | Constr_regexp r -> add r.rexp_trm_kind
    | _ -> ()
    end;
    ignore (constr_map explore c);
    c
    in
  let iter_in_target tg =
    List.iter (fun c -> ignore (explore c)) tg in
  List.iter iter_in_target tgs;
  !res

(* [get_target_regexp_kinds tgs: gets the list of the regexp characterizing
   toplevel functions that appear in the targets that contain constraints based
   on string representation. If one of the targets does not contain the subsequence
   [cTop name] or [cTopFun name], which generate
   [Constr_target [Constr_root; Constr_depth (DepthAt 1); Constr_decl_fun (Some rexp)]],
   then the result will be [None]. *)

exception Topfuns_cannot_filter
let get_target_regexp_topfuns_opt (tgs : target list) : constr_name list option =

  let has_regexp (c : constr) : bool =
    let answer = ref false in
      let rec aux c = (* LATER: optimize using a constr_iter instead of constr_map *)
        match c with
        | Constr_regexp _ -> answer := true; c
        | _ -> ignore (constr_map aux c); c
        in
      ignore (aux c);
      !answer in
  (* Printf.printf "get_target_regexp_topfuns_opt %d\n" (List.length tgs); *)
  let tgs = List.filter (fun tg -> List.exists has_regexp tg) tgs in
  (*Printf.printf "get_target_regexp_topfuns_opt filter %d\n" (List.length tgs);*)
  try
    let constr_names : constr_name list ref = ref [] in
    let rec find_in_target (cs : constr list) : unit =
      match cs with
      | Constr_target [ Constr_root;
                        Constr_depth (DepthAt 1);
                        Constr_decl_fun (_, ((Some _) as constr_name), _, _,_) ]
            :: _ ->
          constr_names := constr_name :: !constr_names
      | _ :: cs2 -> find_in_target cs2
      | [] -> raise Topfuns_cannot_filter
      in
    List.iter find_in_target tgs;
    (*Printf.printf "get_target_regexp_topfuns_opt Some %d\n" (List.length !constr_names);*)
    Some !constr_names
  with Topfuns_cannot_filter ->
    (* Printf.printf "get_target_regexp_topfuns_opt None\n";*)
    None


(******************************************************************************)
(*                        Preprocessing of targets before resolution          *)
(******************************************************************************)

(* [target_flatten tg]: flattens all the constraints of type Constr_target *)
let target_flatten (tg : target) : target =
    let rec aux (cs : target) : target =
      match cs with
      | [] -> []
      | c::cs2 ->
          let r = match c with
            | Constr_target cs1 -> (aux cs1)
            | _ -> [c]
            in
          r @ (aux cs2)
      in
    aux tg

(* [target_to_target_struct]: converts a target into a target struct  *)
let target_to_target_struct (tr : target) : target_struct =
  let tr = target_flatten tr in
  let relative = ref None in
  let occurences = ref None in
  let process_constr (c : constr) : unit =
    match c with
    | Constr_relative re ->
      begin match !relative with
      | None -> relative := Some re;
      | Some _ -> fail None  "Constr_relative provided twice in path"
      end
    | Constr_occurrences oc ->
      begin match !occurences with
      | None -> occurences := Some oc;
      | _ -> fail None "Constr.Constr_occurrences provided twice in path"
      end
    | _ -> ()
    in
  List.iter process_constr tr;
  let tgs = {
    target_path = List.filter (function | Constr_relative _ | Constr_occurrences _ -> false | _ -> true) tr;
    target_relative = begin match !relative with | None -> TargetAt | Some re -> re end;
    target_occurrences = begin match !occurences with | None -> ExpectedOne | Some oc -> oc end; } in
  tgs
(* [is_target_between tr]: checks if [tr] contains a relative constraint different from TargetAt *)
let is_target_between (tr : target) : bool =
   let tgs = target_to_target_struct tr in
   tgs.target_relative <> TargetAt

(******************************************************************************)
(*                              Target resolution                             *)
(******************************************************************************)

(*
  Particular case for target resolution: heap allocated variables
  Patterns:
    - declaration:A declaration is heap allocated if it is mutable
    - usage: particular use of get/access because variables are pointers
      in practice, only dereferencing has to be taken into account: array and
      struct get/access are better matched with regexp
  The user expresses targets as if the variables were not heap allocated but target
  resolution computes the appropriate target
 *)

(*
  Other particular case: accesses
  Pattern:
    when a Rvalue is expected, there is a star operator at the root of the
    subterm
  The user should ignore the existence of this star operator but target resolution
  should compute the appropriate target
 *)

(*
  translate t into the trm the user sees by removing heap allocation
  should be locally applied to the patterns described above
 *)

(* [add_dir d dll]: extends current explicit paths with a direction *)
let add_dir (d : dir) (dll : paths) : paths =
  List.map (fun dl -> d :: dl) dll

(* [is_equal_lit l l']: compares literals l and l' based on their values *)
let is_equal_lit (l : lit) (l' : lit) =
  match l, l' with
  | Lit_unit, Lit_unit -> true
  | Lit_uninitialized, Lit_uninitialized -> true
  | Lit_bool b, Lit_bool b' when b = b' -> true
  | Lit_int n, Lit_int n' when n = n' -> true
  | Lit_double d, Lit_double d' when d = d' -> true
  | Lit_string s, Lit_string s' when s = s' -> true
  | Lit_nullptr, Lit_nullptr -> true
  | _ -> false

(* [get_trm_kind t]: gets the kind of trm [t] *)
(* LATER: we may want to save the kind inside the term? *)
let get_trm_kind (t : trm) : trm_kind =
   let is_unit = begin match t.typ with
                 | Some ty ->
                    begin match ty.typ_desc with
                    | Typ_unit -> true
                    | _ -> false
                    end
                 | None -> false
                 end
    in
   match t.desc with
   | Trm_val _ -> if is_unit then TrmKind_Instr else TrmKind_Expr
   | Trm_var _ -> TrmKind_Expr
   | Trm_record _ | Trm_array _ -> TrmKind_Expr
   | Trm_let_fun _ -> TrmKind_Ctrl (* purposely not an instruction *)
   | Trm_let _ | Trm_let_mult _ -> TrmKind_Instr
   | Trm_typedef _ | Trm_let_record _-> TrmKind_Typedef
   | Trm_if _-> if is_unit then TrmKind_Ctrl else TrmKind_Expr
   | Trm_seq _ -> TrmKind_Ctrl
   | Trm_apps _ -> if is_unit then TrmKind_Instr else TrmKind_Expr
   | Trm_while _ | Trm_do_while _ | Trm_for_c _ | Trm_for _| Trm_switch _ | Trm_abort _ | Trm_goto _ -> TrmKind_Ctrl
   | Trm_omp_routine _ | Trm_extern _  | Trm_namespace _ | Trm_template _ | Trm_arbitrary _ | Trm_using_directive _ -> TrmKind_Any


(* [match_regexp_str r s]: checks if [s] can be matched with [r] *)
let match_regexp_str (r : rexp) (s : string) : bool =
  (*if s = "x" then incr Debug.counter;
  if !Debug.counter = 2 then raise Debug.Breakpoint; *)
  if r.rexp_substr then begin
    try let _ = Str.search_forward r.rexp_exp s 0 in true
    with Not_found -> false
  end else begin
    (* Here we assume that [r.rexp_exp] ends with [^], to ensure that
       the pattern matches all of [s] and not a strict substring of [s]. *)
    Str.string_match r.rexp_exp s 0
  end


(* LATER: this could be passed as argument throughout the function calls *)
(* This variable should only be modified by [Target.with_stringreprs_available] *)
(* Keep in mind that the stringrepr is not available for AST nodes that are removed
   during the computation of [cfeatures_intro], that is, during the translations
   from OptiTrust AST to the C AST. *)
let stringreprs : (stringreprid, string) Hashtbl.t option ref = ref None

(* [print_stringreprs ()]: for debugging purpose *)
let print_stringreprs () : unit =
  match !stringreprs with
  | None -> fail None "Constr.print_stringreprs: no table registered"
  | Some m ->
      let pr id s =
        Printf.printf "stringreprs[%d] = %s\n----\n" id s in
      Printf.printf "====<constr.stringreprs>====\n";
      Hashtbl.iter pr m;
      Printf.printf "====</constr.stringreprs>====\n"

(* [get_stringrepr t]: returns the string representation saved in table [stringreprs],
   or an empty string otherwise *)
let get_stringrepr (t : trm) : string =
    match !stringreprs with
    | None -> fail t.loc (Printf.sprintf "Constr.get_stringrepr: stringreprs must be computed and registered before resolving constraints, %s" (Ast_to_text.ast_to_string t))
    | Some m ->
        match Ast.trm_get_stringreprid t with
        | Some id ->
          begin match Hashtbl.find_opt m id with
          | None -> ""
              (* This term must correspond to a node that was removed during
                 [cfeatures_intro], hence not printed *)
              (* FOR debug: print_stringreprs(); AstC_to_c.trm_print_debug t; *)
          | Some s -> s
          end
        | None -> ""

(* [match_regexp_trm_kind k t]: checks if [t] is of kind [k] *)
let match_regexp_trm_kind (k : trm_kind) (t : trm) : bool =
  (k = TrmKind_Any) || (get_trm_kind t = k)

(* [match_regexp_trm_kinds ks t]: checks if [t] is of one of the kinds in [ks] *)
let match_regexp_trm_kinds (ks : trm_kind list) (t : trm) : bool =
  List.mem (get_trm_kind t) ks

(* [match_regexp_trm r t]: checks if the string representation of [t] can be matched with [r] *)
let match_regexp_trm (r : rexp) (t : trm) : bool =
  if not (match_regexp_trm_kind r.rexp_trm_kind t) then false else begin
    let s = get_stringrepr t in
    (* FOR DEBUG: Printf.printf "Considered: %s\n" s; *)
    s <> "" && match_regexp_str r s
    (* If the stringrepr is not available, we return false *)
  end

(* [is_constr_regexp c]: checks if [c] is a regexp constraint *)
let is_constr_regexp (c : constr) : bool =
  match c with | Constr_regexp _ -> true | _ -> false


(* [check_hastype pred t]: tests whether [t] carries a type
   that satisfies [pred]. If [t] does not carry a type information,
   then the return value is [false]. For constraints to work properly,
   one may want to check that types are up to date (e.g. by reparsing). *)
let check_hastype (pred : typ->bool) (t : trm) : bool =
    match t.typ with
    | Some ty -> pred ty
    | None -> false

(* [check_constraint c]: checks if constraint c is satisfied by trm t *)
let rec check_constraint (c : constr) (t : trm) : bool =
   if trm_has_cstyle Multi_decl t then
     (*
       check the constraint on each element of the seq and return true if one
       is true
      *)
     begin match t.desc with
     | Trm_seq tl -> List.mem true (List.map (check_constraint c) (Mlist.to_list tl))
     | _ -> fail t.loc "Constr.check_constraint: bad multi_decl annotation"
     end
  else

     let _loc = t.loc in
     begin match c, t.desc with
     (*
       target constraints never hold since they are checked against nodes before
       calling check_constraint in resolve_target
      *)
      | Constr_depth _,_
       | Constr_dir _, _
       | Constr_include _, _ ->
        false
     | Constr_regexp r, _ -> match_regexp_trm r t
     | Constr_for_c (p_init, p_cond, p_step, p_body),
       Trm_for_c (init, cond, step, body) ->
        check_target p_init init &&
        check_target p_cond cond &&
        check_target p_step step &&
        check_target p_body body
     | Constr_for (p_index, p_start, p_direction, p_stop, p_step, p_body), Trm_for(l_range, body) ->
        let (index, start, direction, stop, step, _is_parallel) = l_range in
        let direction_match = match p_direction with
        | None -> true
        | Some d -> d = direction in
        check_name p_index index &&
        direction_match &&
        check_target p_start start &&
        check_target p_stop stop &&
        check_target p_step (loop_step_to_trm step) &&
        check_target p_body body
     | Constr_while (p_cond, p_body), Trm_while (cond, body) ->
        check_target p_cond cond &&
        check_target p_body body
     | Constr_do_while (p_body, p_cond), Trm_do_while (body, cond) ->
        check_target p_body body &&
        check_target p_cond cond
     | Constr_if (p_cond, p_then, p_else), Trm_if (cond, then_t, else_t) ->
        check_target p_cond cond &&
        check_target p_then then_t &&
        check_target p_else else_t
      | Constr_decl_var (ty_pred, name, p_body) , Trm_let (_,(x,tx), body) ->
        ty_pred (get_inner_ptr_type tx) &&
        check_name name x &&
        check_target p_body body
     | Constr_decl_fun (ty_pred, name, cl_args, p_body,is_def),
       Trm_let_fun (x, tx, args, body) ->
        let body_check =
          let is_body_unit = is_trm_uninitialized body in
          if is_def then (check_target p_body body && not (is_body_unit))
           else is_body_unit in
        ty_pred tx &&
        check_name name x &&
        check_args cl_args args &&
        body_check
     | Constr_decl_type name, Trm_typedef td ->
        let is_new_typ = begin match td.typdef_body with
        | Typdef_alias _ -> true
        | Typdef_record _ -> true
        | _ -> false
        end in
        let x = td.typdef_tconstr in
        is_new_typ && check_name name x
     | Constr_decl_enum (name, cec), Trm_typedef td ->
        begin match td.typdef_body with
        | Typdef_enum xto_l -> check_name name td.typdef_tconstr && check_enum_const cec xto_l
        | _ -> false
        end
     | Constr_seq cl, Trm_seq tl when
        not ((List.exists (function No_braces _ -> true | _ -> false) t.annot.trm_annot_cstyle) || List.mem Main_file t.annot.trm_annot_files)->
        check_list ~depth:(DepthAt 0) cl (Mlist.to_list tl) (* LATER/ check why depth 0 here and not
        in constra_app *)
     | Constr_var name, Trm_var (_, x) ->
        check_name name x
     | Constr_lit pred_l, Trm_val (Val_lit l) ->
        pred_l l
     | Constr_app (p_fun, cl_args, accept_encoded), Trm_apps (f, args) ->
        if not accept_encoded then
          begin match f.desc with
          | Trm_val (Val_prim (Prim_new _))
          | Trm_val (Val_prim (Prim_unop Unop_get)) -> false
          |  _ -> check_target p_fun f &&
                  check_list ~depth:(DepthAny) cl_args args
          end
        else
          check_target p_fun f &&
          check_list ~depth:(DepthAny) cl_args args
     | Constr_label (so, p_body), _ ->
        let t_labels = trm_get_labels t in
        List.fold_left (fun acc l -> check_name so l || acc) false t_labels &&
        check_target p_body t
     | Constr_goto so, Trm_goto l ->
        check_name so l
     | Constr_return p_res, Trm_abort (Ret (Some res)) ->
        check_target p_res res
     | Constr_abort Any, Trm_abort _ -> true
     | Constr_abort Return, Trm_abort (Ret _) -> true
     | Constr_abort Break, Trm_abort (Break _) -> true
     | Constr_abort Continue, Trm_abort (Continue _) -> true
     | Constr_access (p_base, ca,ia), _ ->
        let (base, al) = get_nested_accesses t in
        check_target p_base base &&
        check_accesses ~inner_accesses:ia ca al
     | Constr_switch (p_cond, cc), Trm_switch (cond, cases) ->
        check_target p_cond cond &&
        check_cases cc cases
     | Constr_bool b, _ -> b
     | Constr_root, _ -> trm_is_mainfile t

     | Constr_prim pred, Trm_val (Val_prim p1) ->
        pred p1
     | Constr_mark (pred, _m), _ ->
        let t_marks = trm_get_marks t in
        begin match t.desc with
        | Trm_seq tl | Trm_array tl ->
          (List.exists pred t_marks) || (List.fold_left (fun acc x -> (List.exists pred x) || acc) false tl.marks)
        | Trm_record tl -> (List.exists pred t_marks) || (List.fold_left (fun acc x -> (List.exists pred x) || acc) false tl.marks)
        | _ -> List.exists pred t_marks
        end
     | Constr_hastype pred , _ ->
        check_hastype pred t
     | Constr_var_init , Trm_apps ({desc = Trm_val (Val_prim (Prim_new _)); _}, [arg]) -> false
     | Constr_var_init, Trm_val (Val_lit (Lit_uninitialized)) -> false
     | Constr_var_init , _ -> true
     | Constr_array_init, Trm_array _ -> true
     | Constr_struct_init, Trm_record _ -> true
     | Constr_omp (pred, _), _ -> trm_has_pragma pred t
     | _ -> false
     end

(* [check_name name s]: checks if constraint [name] matches string [s] *)
and check_name (name : constr_name) (s : string) : bool =
  match name with
  | None -> true
  | Some r ->
     match_regexp_str r s

(* [check_list ~depth lpred tl]: checks if [tl] satisfy the predicate [lpred] *)
and check_list ?(depth : depth = DepthAny) (lpred : target_list_pred) (tl : trms) : bool =
  let ith_target = lpred.target_list_pred_ith_target in
  let validate = lpred.target_list_pred_validate in
  validate (List.mapi (fun i t -> check_target ~depth (ith_target i) t) tl)

(* [check_args lpred tx]: checks if [txl] satisfy the predicate [lpred] *)
and check_args (lpred : target_list_pred) (txl : typed_vars) : bool =
  let ith_target = lpred.target_list_pred_ith_target in
  let validate = lpred.target_list_pred_validate in
  validate (List.mapi (fun i tx -> check_arg (ith_target i) tx) txl)

(* [check_arg tg tx]: checks if the typed-variable [tx] satisfies the argument-constraint [tg].
   An argument constraint is a singleton constraint, of the form [cArg ..] or [cTrue]. *)
and check_arg (tg:arg_constraint) ((var_name, var_typ) : typed_var) : bool =
  match tg with
  | [] -> true
  | [c] -> begin match c with
           | Constr_bool true -> true
           | Constr_arg (var_constraint, typ_constraint) ->
              var_constraint var_name && typ_constraint var_typ
           | _ -> fail None "Constr.check_arg: target expressing constraints on arguments must be of
                             the form [cArg...] or [cTrue]."
          end
  | _ -> fail None "Constr.check_arg: target expressing constraints on arguments must be list with at most one item."


(* [check_accesses ~inner_accesses ca al]: checks if [al] is matched by [ca] *)
and check_accesses ?(inner_accesses : bool = true) (ca : constr_accesses) (al : trm_access list) : bool =
  let rec aux (cal : constr_access list) (al : trm_access list) : bool =
    match cal, al with
    | [], a -> if not inner_accesses
                  then begin match a with
                       | [] -> true
                       | _  -> false
                       end
                  else true
    | Array_access p_index :: cal, Array_access_get index :: al ->
       check_target p_index index &&
       aux cal al
    | Array_access p_index :: cal, Array_access_addr index :: al ->
       check_target p_index index &&
       aux cal al
    | Struct_access so :: cal, Struct_access_get f :: al ->
       check_name so f &&
       aux cal al
    | Struct_access so :: cal, Struct_access_addr f :: al ->
       check_name so f &&
       aux cal al
    | Any_access :: cal, _ :: al -> aux cal al
    | _ -> false
  in
  match ca with
  | None -> true
  | Some cal -> aux cal al

(* [check_cases cc cases]: checks if [cases] are matched by [cc] *)
and check_cases (cc : constr_cases) (cases : (trms * trm) list) : bool =
  let rec aux (cl : (case_kind * target) list)
    (cases : (trms * trm) list) : bool =
    match cl, cases with
    | [], [] -> true
    | (k, p_body) :: cl, (tl, body) :: cases ->
       check_kind k tl &&
       check_target p_body body &&
       aux cl cases
    | _ -> false
  in
  match cc with
  | None -> true
  | Some cl -> aux cl cases

(* [check_kind k tl]: checks if [tl] are of kind [k] *)
and check_kind (k : case_kind) (tl : trms) : bool =
  match k, tl with
  | Case_any, _ -> true
  | Case_default, [] -> true
  | Case_val p_val, _ ->
     (* check that one of the cases corresponds to the given target *)
     List.mem true (List.map (check_target p_val) tl)
  | _ -> false

(* [check_enum_const cec xto_l]: checks if [xto_l] are matched by constraint [cec] *)
and check_enum_const (cec : constr_enum_const)
  (xto_l : (string * (trm option)) list) : bool =
  match cec with
  | None -> true
  | Some cnp_l ->
     List.for_all2
       (fun (cn, p) (n, t_o) ->
         check_name cn n &&
         match p, t_o with
         | [], None -> true
         | _, None -> false
         | _, Some t -> check_target p t
       )
       cnp_l
       xto_l

(* [check_target ~depth tr t]: checks if target tr leads to at least one subterm of t *)
and check_target ?(depth : depth = DepthAny) (tr : target) (t : trm) : bool =
  match resolve_target_simple ~depth tr t with
  | [] -> false
  | _ -> true

(*
  resolve_target computes the directions to matching subterms
  expected invariants: no duplicate target and no target which is prefix of
  another target that appears after it in the list. Guaranteed by the call to
  sort_unique *)

(* [debug_resolution]: only for debugging *)
and debug_resolution = false

(* [resolve_target_simple ~depth trs t]: the main function that resolves a target
     It returns a list of paths that resolve to the given target
     [depth] is the depth level where the targetd should be resolved
     [trs]: a clean target t
     [t]: ast  *)
and resolve_target_simple ?(depth : depth = DepthAny) (trs : target_simple) (t : trm) : paths =
  Stats.incr_target_resolution_steps ();
  let epl =
    match trs with
    | [] -> [[]]
    | Constr_target tl :: trest -> (* LATER: see if we can flatten recursively in targets before we start *)
       resolve_target_simple ~depth (tl @ trest) t

    | Constr_or tl :: trest -> (* LATER: maybe we'll add an option to enforce that each target from the list tl resolves to at least one solution *)
        let all_targets_must_resolve = false in
        let res = List.fold_left (fun acc tr ->
          let potential_targets = resolve_target_simple ~depth (tr @ trest) t in
          begin match potential_targets with
          | ([] | [[]]) when all_targets_must_resolve -> fail t.loc "Constr.resolve_target_simple: for Constr_and all targets should match a trm"
          | _ ->
            Path.union acc potential_targets
          end ) [] tl in
       if debug_resolution then begin
          Printf.printf "resolve_target_simple[Constr_or]\n  ~target:%s\n  ~term:%s\n  ~res:%s\n"
            (target_to_string trs)
            (AstC_to_c.ast_to_string  t)
            (paths_to_string ~sep:"\n   " res)
        end;
        res
    | Constr_diff (tl1 , tl2) :: [] ->
      let all_targets_must_resolve = false in
      let targets_to_keep =
      List.fold_left (fun acc tr ->
          let potential_targets = resolve_target_simple ~depth tr t in
          begin match potential_targets with
          | ([] | [[]]) when all_targets_must_resolve -> fail t.loc "Constr.resolve_target_simple: for Constr_and all targets should match a trm"
          | _ ->
            Path.union acc potential_targets
          end ) [] tl1 in
      let targets_to_remove =
      List.fold_left (fun acc tr ->
          let potential_targets = resolve_target_simple tr t in
          begin match potential_targets with
          | ([] | [[]]) -> acc
          | _ ->
            Path.union acc potential_targets
          end ) [] tl2 in
      Path.diff targets_to_keep targets_to_remove

    | Constr_and tl :: trest ->
        (* LATER: ARTHUR : optimize resolution by resolving the targets only by exploring
          through the paths that are candidates; using e.g. path_satisfies_target *)
        let all_targets_must_resolve = false in
        Xlist.fold_lefti (fun i acc tr ->
          let targetsi = resolve_target_simple (tr @ trest) t in
          begin match targetsi with
          | ([] | [[]]) when all_targets_must_resolve -> fail t.loc "Constr.resolve_target_simple: for Constr_and all targets should match a trm"
          | _ ->
            if i = 0
              (* First step, initalize the acc *)
              then targetsi
            (* Compute the intersection of all resolved targets *)
              else Path.intersect acc targetsi
          end) [] tl

    | Constr_diff (_ , _) :: _ ->
        fail t.loc "Constr.cDiff can only appear at last item in a target"

    | Constr_depth new_depth :: tr ->
        (* Force the depth argument for the rest of the target, override the current [depth] *)
        resolve_target_simple ~depth:new_depth tr t
    | Constr_dir d :: tr ->
        follow_dir d tr t
    | c :: p ->
      let strict = match depth with
        | DepthAt 0 -> true
        | _ when c = Constr_root -> true (* useless to search for the root in depth *)
        | _ -> false in
      let skip_here = match depth with DepthAt n when n > 0 -> true | _ -> false in
      let res_deep =
        if strict
           then [] (* in strict mode, must match c here *)
           else (explore_in_depth ~depth (c :: p) t) in
      let res_here =
         if skip_here || (is_constr_regexp c && res_deep <> [])
           then [] (* if a regexp matches in depth, don't test it here *)
           else (resolve_constraint c p t) in

      (* DEBUG *)
      if debug_resolution then begin
         Printf.printf "resolve_target_simple\n  ~strict:%s\n  ~target:%s\n  ~term:%s\n ~deep:%s\n  ~here:%s\n"
          (if strict then "true" else "false")
          (target_to_string trs)
          (AstC_to_c.ast_to_string  t)
          (paths_to_string ~sep:"\n   " res_deep)
          (paths_to_string ~sep:"\n   " res_here);
      end;
        (* Printf.printf " ~deep:%s\n  ~here:%s\n"
          (paths_to_string ~sep:"\n   " res_deep)
          (paths_to_string ~sep:"\n   " res_here); *)

      res_deep@res_here  (* put deeper nodes first *) in
  List.sort_uniq compare_path epl

(* [resolve_target_struct tgs t]: resolves the structured target [tgs] over trm [t]
    This function gets the result from [resolve_target_simple] and checks if it matches
    the given requirements. If one of the requirements is not fulfilled the target resolution will fail. *)
and resolve_target_struct (tgs : target_struct) (t : trm) : paths =
  let res = resolve_target_simple tgs.target_path t in
  let nb = List.length res in
  (* Check if nb is equal to the specification of tgs.target_occurrences, if not then something went wrong *)
  let error s =
    raise (Resolve_target_failure (None, s)) in
  begin match tgs.target_occurrences with
  | ExpectedOne -> if nb <> 1 then error (Printf.sprintf "resolve_target_struct: expected exactly one match, got %d." nb); res
  | ExpectedNb n -> if nb <> n then error (Printf.sprintf "resolve_target_struct: expected %d matches, got %d." n nb); res
  | ExpectedMulti -> if nb = 0 then error (Printf.sprintf "resolve_target_struct: expected at least one occurrence, got %d." nb); res
  | ExpectedAnyNb -> res
  | ExpectedSelected (n_opt, i_selected) ->
    begin match n_opt with
    | Some n ->
      if n <> nb then error (Printf.sprintf "resolve_target_struct: expected %d matches, got %d" n nb)
        else
          Xlist.filter_selected i_selected res;
    | None -> if nb = 0
                then error (Printf.sprintf "resolve_target_struct: expected %d matches, got %d" (List.length i_selected) nb)
                else Xlist.filter_selected i_selected res;
    end
  | FirstOcc -> [fst (Xlist.uncons res)]
  | LastOcc ->  [snd (Xlist.unlast res)]
  end

(* [resolve_target tg]: resolves the target [tg] *)
and resolve_target (tg : target) (t : trm) : paths =
  let tgs = target_to_target_struct tg in
  if tgs.target_relative <> TargetAt
    then fail None "Constr.resolve_target: this target should not contain a tBefore/tAfter/tFirst/tLast";
  try resolve_target_struct tgs t
  with Resolve_target_failure (_loc_opt,str) ->
    fail None ("Constr." ^ str ^ "\n" ^ (target_to_string tg))

(* [resolve_target_exactly_one tg t]: similar to [resolve_target] but this one fails if the target resolves to more than one path *)
and resolve_target_exactly_one (tg : target) (t : trm) : path =
  match resolve_target tg t with
  | [p] -> p
  | _ -> fail t.loc "Constr.resolve_target_exactly_one: obtained several targets."

(* [resolve_constraint c p t]: checks [c] against [t] and in case of success continue with [p] *)
and resolve_constraint (c : constr) (p : target_simple) (t : trm) : paths =
  let loc = t.loc in
  match c with
  | Constr_target _ -> fail t.loc "Constr.resolve_constraint should not reach a Constr_target"
  (*
    do not resolve in included files, except if the constraint is Constr_include
   *)
  | Constr_include h when trm_is_include t ->
     (* remove the include annotation for target resolution to proceed in the
       included file *)
     resolve_target_simple p (trm_alter ~annot:(Some trm_annot_default) t)
  | _ when trm_is_include t ->
     print_info loc "Constr.resolve_constraint: not an include constraint\n";
     []
  (* target constraints first *)
  (* following directions *)
  | Constr_dir d -> follow_dir d p t
  (*
    if the constraint is a target constraint that does not match the node or
    if it is another kind of constraint, then we check if it holds
   *)
  | c when check_constraint c t -> resolve_target_simple p t
  | _ ->
     print_info loc "Constr.resolve_constraint: constraint %s does not hold\n"
       (constr_to_string c);
     []

(* [explore_in_depth ~depth p t]: calls resolve_target_simple on subterms of t if possible *)
and explore_in_depth ?(depth : depth = DepthAny) (p : target_simple) (t : trm) : paths =
  (* By default, traversing a node decreases the depth *)
  let aux = resolve_target_simple ~depth:(depth_pred depth) p in
  (* For bodies of functions/loops/conditials, however, we do decrease the depth
     because traversing the [Trm_seq] just below will already decrease the depth *)
  let aux_body = resolve_target_simple ~depth p in

  let loc = t.loc in
  (* no exploration in depth in included files *)
  if trm_is_include t then begin
     print_info loc "Constr.explore_in_depth: no exploration in included files\n";
     []
     end
  else if trm_has_cstyle Multi_decl t then
     (* explore each declaration in the seq *)
     begin match t.desc with
     | Trm_seq tl ->
        explore_list (Mlist.to_list tl) (fun i -> Dir_seq_nth i) (explore_in_depth p)
     | _ -> fail loc "Constr.explore_in_depth: bad multi_decl annotation"
     end
  else
     begin match t.desc with
     | Trm_let (_ ,(_, _), body) ->
       add_dir Dir_body (aux body)
     | Trm_let_fun (_, _ , _,body) ->
        add_dir Dir_body (aux_body body)
     | Trm_typedef td  ->
      begin match td.typdef_body with
      | Typdef_enum xto_l ->
        let (il, tl) =
          Xlist.fold_lefti
            (fun n (il, tl) (_, t_o) ->
              match t_o with
              | None -> (il, tl)
              | Some t -> (il@[n], tl@[t])
            )
           ([], [])
           xto_l
        in
        add_dir Dir_name (aux (trm_var ~loc td.typdef_tconstr)) @
        (explore_list (List.map (fun (y, _) -> trm_var ~loc y) xto_l)
           (fun n -> Dir_enum_const (n, Enum_const_name))
           (aux)) @
        (explore_list tl
           (fun n -> Dir_enum_const (List.nth il n, Enum_const_val))
           (aux))
      | _ -> []
      end
     | Trm_abort (Ret (Some body)) ->
        add_dir Dir_body (aux body)
     | Trm_for (l_range, body) ->
        let ( _, start, _, stop, step, is_parallel) = l_range in
        let step_t = loop_step_to_trm step in
        (add_dir Dir_for_start (aux start)) @
        (add_dir Dir_for_stop (aux stop)) @
        (add_dir Dir_for_step (aux step_t)) @
        (add_dir Dir_body (aux_body body))
     | Trm_for_c (init, cond, step, body) ->
        (* init *)
        (add_dir Dir_for_c_init (aux init)) @
        (* cond *)
        (add_dir Dir_cond (aux cond)) @
        (* step *)
        (add_dir Dir_for_c_step (aux step)) @
        (* body *)
        (add_dir Dir_body (aux_body body))
     | Trm_while (cond, body) ->
        (* cond *)
        (add_dir Dir_cond (aux cond)) @
        (* body *)
        (add_dir Dir_body (aux_body body))
     | Trm_do_while (body, cond) ->
        (* body *)
        (add_dir Dir_body (aux_body body)) @
        (* cond *)
        (add_dir Dir_cond (aux cond))
     | Trm_if (cond, then_t, else_t) ->
        (* cond *)
        (add_dir Dir_cond (aux cond)) @
        (* then *)
        (add_dir Dir_then (aux_body then_t)) @
        (* else *)
        (add_dir Dir_else (aux_body else_t))
     | Trm_apps (f, args) ->
        (* fun *)
        (add_dir Dir_app_fun (aux f)) @
        (* args *)
        (explore_list args (fun n -> Dir_arg_nth n) (aux))
     | Trm_array tl ->
        explore_list (Mlist.to_list tl) (fun n -> Dir_array_nth n) (aux)
     | Trm_seq tl ->
        explore_list (Mlist.to_list tl) (fun n -> Dir_seq_nth n) (aux)
     | Trm_record tl ->
        explore_list (Xlist.split_pairs_snd (Mlist.to_list tl)) (fun n -> Dir_struct_nth n) (aux)
     | Trm_switch (cond, cases) ->
        (add_dir Dir_cond (aux cond)) @
        (Xlist.fold_lefti (fun i epl case -> epl@explore_case depth i case p) [] cases)
     | _ ->
        print_info loc "explore_in_depth: cannot find a subterm to explore\n";
        []
     end

(*
  call resolve_target_simple on given case name and body
  i is the index of the case in its switch
 *)
(* [explore_case depth i case p]: calls resolve_target_simple on given case name and body,
    [i] is the index of the case in its switch satement *)
and explore_case (depth : depth) (i : int) (case : trms * trm) (p : target_simple) : paths =
  let aux = resolve_target_simple ~depth p in
  let (tl, body) = case in
  match tl with
  (* default case *)
  | [] ->
     add_dir (Dir_case (i, Case_body)) (aux body)
  | _ ->
     (Xlist.fold_lefti
        (fun j epl t ->
          epl @
          (add_dir (Dir_case (i, Case_name j)) (aux t))
        )
        []
        tl
     ) @
     add_dir (Dir_case (i, Case_body)) (aux body)

(* [follow_dir d p t]: follows the direction [d] in [t] and call resolve_target_simple on [p] *)
and follow_dir (d : dir) (p : target_simple) (t : trm) : paths =
  let aux = resolve_target_simple p in
  let loc = t.loc in
  match d, t.desc with
  | Dir_array_nth n, Trm_array tl ->
    app_to_nth_dflt loc (Mlist.to_list tl) n
       (fun nth_t -> add_dir (Dir_array_nth n) (aux nth_t))
  | Dir_seq_nth n, Trm_seq tl ->
    app_to_nth_dflt loc (Mlist.to_list tl) n
       (fun nth_t -> add_dir (Dir_seq_nth n) (aux nth_t))
  | Dir_struct_nth n, Trm_record tl ->
     app_to_nth_dflt loc (Xlist.split_pairs_snd (Mlist.to_list tl)) n
       (fun nth_t -> add_dir (Dir_struct_nth n) (aux nth_t))
  | Dir_cond, Trm_if (cond, _, _)
    | Dir_cond, Trm_while (cond, _)
    | Dir_cond, Trm_do_while (_, cond)
    | Dir_cond, Trm_for_c (_, cond, _, _)
    | Dir_cond, Trm_switch (cond, _) ->
     add_dir Dir_cond (aux cond)
  | Dir_then, Trm_if (_, then_t, _) ->
     add_dir Dir_then (aux then_t)
  | Dir_else, Trm_if (_, _, else_t) ->
     add_dir Dir_else (aux else_t)
  | Dir_var_body, Trm_let (_, _, body) ->
     let new_op_arg = new_operation_arg body in
     add_dir Dir_var_body (aux new_op_arg)
  | Dir_body, Trm_let (_, _,body)
    | Dir_body, Trm_let_fun (_, _, _, body)
    | Dir_body, Trm_for_c (_, _, _, body)
    | Dir_body, Trm_for (_, body)
    | Dir_body, Trm_while (_, body)
    | Dir_body, Trm_do_while (body, _)
    | Dir_body, Trm_abort (Ret (Some body)) ->
     add_dir Dir_body (aux body)
  | Dir_for_c_init, Trm_for_c (init, _, _, _) ->
     add_dir Dir_for_c_init (aux init)
  | Dir_for_c_step, Trm_for_c (_, _, step, _) ->
     add_dir Dir_for_c_step (aux step)
  | Dir_for_start, Trm_for (l_range, _) ->
     let (_, start,  _, _, _, _) = l_range in
     add_dir Dir_for_start (aux start)
  | Dir_for_stop, Trm_for (l_range, _) ->
     let (_, _, _, stop, _, _) = l_range in
     add_dir Dir_for_stop (aux stop)
  | Dir_app_fun, Trm_apps (f, _) -> add_dir Dir_app_fun (aux f)
  | Dir_arg_nth n, Trm_apps (_, tl) ->
     app_to_nth_dflt loc tl n (fun nth_t ->
         add_dir (Dir_arg_nth n) (aux nth_t))
  | Dir_arg_nth n, Trm_let_fun (_, _, arg, _) ->
     let tl = List.map (fun (x, _) -> trm_var ~loc x) arg in
     app_to_nth_dflt loc tl n (fun nth_t ->
         add_dir (Dir_arg_nth n) (aux nth_t))
  | Dir_name, Trm_typedef td ->
     add_dir Dir_name (aux (trm_var ~loc td.typdef_tconstr))
  | Dir_name, Trm_let (_,(x,_),_)
    | Dir_name, Trm_let_fun (x, _, _, _)
    | Dir_name, Trm_goto x ->
     add_dir Dir_name (aux (trm_var ~loc x))
  | Dir_case (n, cd), Trm_switch (_, cases) ->
     app_to_nth_dflt loc cases n
       (fun (tl, body) ->
         match cd with
         | Case_body ->
            add_dir (Dir_case (n, cd)) (aux body)
         | Case_name i ->
            app_to_nth_dflt loc tl i (fun ith_t ->
                add_dir (Dir_case (n, cd)) (aux ith_t))
       )
  | Dir_enum_const (n, ecd), Trm_typedef td ->
     begin match td.typdef_body with
     | Typdef_enum xto_l ->
          app_to_nth_dflt loc xto_l n
          (fun (x, t_o) ->
            match ecd with
            | Enum_const_name ->
               add_dir (Dir_enum_const (n, ecd)) (aux (trm_var ~loc x))
            | Enum_const_val ->
               begin match t_o with
               | None ->
                  print_info loc "follow_dir: no value for constant of index %d\n"
                    n;
                  []
               | Some t ->
                  add_dir (Dir_enum_const (n, ecd)) (aux t)
               end
          )
      | _ -> []
      end
  | _, _ ->
     print_info loc "follow_dir: direction %s does not match"
       (dir_to_string d);
     []

(* [explore_list tl d cont]: calls [cont] on each element of the list and gathers the results
    used for seq, array, struct, and fun arguments
   [d] is the function that gives the direction to add depending on the index *)
and explore_list (tl : trms) (d : int -> dir)
  (cont : trm -> paths) : paths =
  Xlist.fold_lefti (fun i epl t -> epl@add_dir (d i) (cont t)) [] tl

(*
  call cont on each element of the list whose index is in the domain and
  gather the results
  d: function that gives the direction to add depending on the index
 *)

(* [explore_list_ind tl d dom cont]: calls [cont] on each element of the list,
   [index] is in the domain and gather results
   [d] is a fucntion that gives direction to add depending on the [index] *)
and explore_list_ind (tl : trms) (d : int -> dir) (dom : int list)
  (cont : trm -> paths) : paths =
  Xlist.fold_lefti
    (fun i epl t ->
      if List.mem i dom then epl@add_dir (d i) (cont t) else epl)
    []
    tl


(******************************************************************************)
(*                          Target-between resolution                         *)
(******************************************************************************)
(* [extract_last_path_item p]: extracts the last direction from a nonempty path *)
let extract_last_path_item (p : path) : dir * path =
  match List.rev p with
  | [] -> raise Not_found
  | d :: p' -> (d, List.rev p')

(* [extract_last_dir p]: extracts the last direction on a sequence *)
let extract_last_dir (p : path) : path * int =
  match List.rev p with
  | [] -> raise Not_found
  | d :: p' ->
    begin match d with
    | Dir_seq_nth i -> (List.rev p', i)
    | _ -> fail None "Constr.extract_last_dir: expected a directory in a sequence"
    end

(* [get_sequence_length t]: gets the number of instructions on a sequence *)
let get_sequence_length (t : trm) : int =
  begin match t.desc with
  | Trm_seq tl -> Mlist.length tl
  | _ -> fail t.loc "Constr.get_sequence_length: expected a sequence"
  end

(* [get_arity_of_seq_at]: gets the arity of a sequence at path [p] *)
let get_arity_of_seq_at (p : path) (t : trm) : int =
  let (d,p') =
    try extract_last_path_item p
    with Not_found -> fail None "Constr.get_arity_of_seq_at: expected a nonempty path"
    in
  match d with
  | Dir_seq_nth _ ->
      let seq_trm = Path.resolve_path p' t in
      get_sequence_length seq_trm
  | Dir_then | Dir_else | Dir_body  ->
      let seq_trm = Path.resolve_path p t in
      get_sequence_length seq_trm
  | _ -> fail None "Constr.get_arity_of_seq_at: expected a Dir_seq_nth, Dir_then, Dir_else or Dir_body as last direction"

(* [compute_relative_index rel t p]: computes the relative index for relative targets different from [TargetAt] *)
let compute_relative_index (rel : target_relative) (t : trm) (p : path) : path * int =
  match rel with
  | TargetAt -> fail None "Constr.compute_relative_index: Didn't expect a TargetAt"
  | TargetFirst -> (p, 0)
  | TargetLast -> (p, get_arity_of_seq_at p t)
  | TargetBefore | TargetAfter ->
      let shift =
         match rel with
         | TargetBefore -> 0
         | TargetAfter -> 1
         | _ -> assert false
         in
      let (d,p') =
        try extract_last_path_item p
        with Not_found -> fail None "Constr.compute_relative_index: expected a nonempty path"
        in
      match d with
      | Dir_seq_nth i -> (p', i + shift)
      | _ -> fail None "Constr.compute_relative_index: expected a Dir_seq_nth as last direction"

(* [resolve_target_between tg t]: resolves a target that points before or after an instruction *)
let resolve_target_between (tg : target) (t : trm) : (path * int) list =
  let tgs = target_to_target_struct tg in
  if tgs.target_relative = TargetAt
    then fail None "Constr.resolve_target_between:this target should contain a tBefore, tAfter, tFirst, or tLast";
  let res = resolve_target_struct tgs t in
  List.map (compute_relative_index tgs.target_relative t) res

(* [resolve_target_between_exactly_one tg t]: similar to [resolve_target_between] but this one fails if the target matches
    to multiple paths *)
let resolve_target_between_exactly_one (tg : target) (t : trm) : (path * int) =
  match resolve_target_between tg t with
  | [(p,i)] -> (p,i)
  | _ -> fail t.loc "Constr.resolve_target_between_exactly_one: obtainer several targets"
