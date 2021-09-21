open Ast
open Tools
open Path

include Constr

type path = Path.path
type paths = path list

type case_dir = Path.case_dir

(* NOW INCLUDED
type constr = Constr.constr
type target = constr list
*)

(* Shorthands for Trace functions *)
include Trace


(******************************************************************************)
(*                        Smart constructors for targets                      *)
(******************************************************************************)


(*
  a smart constructor builds a target
  thus, the user provides a target using them
  this list is then flattened to call resolve_target
  unit args are used because of optional arguments
  *)

(******************************************************************************)
(*                             Logic constraints                              *)
(******************************************************************************)

let bTrue : constr =
  Constr_bool true

let bFalse : constr =
  Constr_bool false

let cStrict : constr =
  Constr_depth (DepthAt 1)

let cChain (cstrs : constr list) : constr =
  Constr_target cstrs

(******************************************************************************)
(*                             Relative targets                               *)
(******************************************************************************)

let tBefore : constr =
  Constr_relative TargetBefore

let tAfter : constr =
  Constr_relative TargetAfter

let tFirst : constr =
  Constr_relative TargetFirst

let tLast : constr =
  Constr_relative TargetLast

(******************************************************************************)
(*                            Number of targets                               *)
(******************************************************************************)

let nbMulti : constr =
  Constr_occurrences ExpectedMulti

let nbAny : constr =
    Constr_occurrences ExpectedAnyNb

let nbExact (nb : int) : constr =
    Constr_occurrences (ExpectedNb nb)

let tIndices ?(nb : int = -1) (indices : int list) : constr =
  let expected_nb = match nb with
    | -1 -> None
    | _ -> Some nb in
  Constr_occurrences (ExpectedSelected (expected_nb, indices)  )

let tIndex ?(nb : int = -1) (index : int) : constr =
  tIndices ~nb [index]


(******************************************************************************)
(*                                Directions                                  *)
(******************************************************************************)

let target_of_path (p : path) : target =
  List.map (fun d -> Constr_dir d) p

let dRoot : constr =
    Constr_root

let dNth (n : int) : constr =
    Constr_dir (Dir_seq_nth n)

let dCond : constr =
    Constr_dir Dir_cond

let dThen : constr =
    Constr_dir Dir_then

let dElse : constr =
    Constr_dir Dir_else

let dBody : constr =
    Constr_dir Dir_body

let dInit : constr =
    Constr_dir Dir_for_c_init

let dStep : constr =
    Constr_dir Dir_for_c_step

let dArg (n : int) : constr =
    Constr_dir (Dir_arg n)

let dName : constr =
    Constr_dir Dir_name

let dDirCase (n : int) (cd : case_dir) : constr =
    Constr_dir (Dir_case (n, cd))

let dCaseName (n : int) : case_dir = Case_name n

let dCaseBody : case_dir = Case_body

let dEnumConst (n : int)
  (ecd : enum_const_dir) : constr =
    Constr_dir (Dir_enum_const (n, ecd))

let dEnumConstName : enum_const_dir = Enum_const_name

let dEnumConstVal : enum_const_dir = Enum_const_val


(* [string_to_rexp regexp substr s trmKind]  transforms a string into a regular expression
    used to match ast nodes based on their code representation.
    [string_to_rexp] - denotes a flag to tell if the string entered is a regular epxression or no
    [substr] - denotes a flag to decide if we should target strings whcih contain string [s] or not
    [trmKind] - denotes the kind of the ast note represented in code by string [s]
*)
let string_to_rexp (regexp : bool) (substr : bool) (s : string) (trmKind : trm_kind) : rexp =
    { rexp_desc = (if regexp then "Regexp" else "String") ^ "(\"" ^ s ^ "\")";
      rexp_exp = (if regexp then Str.regexp s else
                       if substr then Str.regexp_string s else Str.regexp ("^" ^ s ^ "$"));
      rexp_substr = substr;
      rexp_trm_kind = trmKind; }

  let string_to_rexp_opt (regexp : bool) (substr : bool) (s : string) (trmKind : trm_kind) : rexp option =
    let res =
      if s = ""
        then None
        else Some (string_to_rexp regexp substr s trmKind)
      in
    res

(******************************************************************************)
(*                             String matching                                *)
(******************************************************************************)

let sInstrOrExpr ?(substr : bool = false) (tk : trm_kind) (s : string) : constr =
  Constr_regexp (string_to_rexp false substr s  tk)

let sInstr ?(substr : bool = true) (s : string) : constr =
  sInstrOrExpr ~substr TrmKind_Instr s

let sExpr ?(substr : bool = true) (s : string)  : constr =
  sInstrOrExpr ~substr TrmKind_Expr s

let sInstrOrExprRegexp (tk : trm_kind) (substr : bool) (s : string) : constr =
  Constr_regexp (string_to_rexp true substr s tk)

let sInstrRegexp ?(substr : bool = false) (s : string) : constr =
  sInstrOrExprRegexp TrmKind_Instr substr s

let sExprRegexp ?(substr : bool = false) (s : string) : constr =
  sInstrOrExprRegexp TrmKind_Expr substr s


(******************************************************************************)
(*                                Ast nodes                                   *)
(******************************************************************************)

let cInclude (s : string) : constr =
    Constr_include s

let cSetVar (x : var) : constr =
  sInstr (x ^ " " ^ "=")

let cArg (n : int) : constr =
  Constr_dir (Dir_arg n)

let cVarDef
  ?(regexp : bool = false) ?(substr : bool = false) ?(body : target = []) (name : string) : constr =
  let ro = string_to_rexp_opt regexp substr name TrmKind_Instr in
  let p_body =  body in
    Constr_decl_var (ro, p_body)

let cFor ?(direction : loop_dir = DirUp) ?(start : target = []) ?(stop : target = []) ?(step : target = []) ?(body : target = []) (index : string) : constr =
  let ro = string_to_rexp_opt false false index TrmKind_Instr in
  Constr_for (ro, direction, start, stop, step, body)

let cFor_c ?(init : target = [])
  ?(cond : target = []) ?(step : target = []) ?(body : target = []) (index : string) : constr =
  let init =
      match index, init with
      | "", [] -> init
      | "", _ -> init
      | _, [] -> [cVarDef index]
      | _, _::_ -> init
      in
    Constr_for_c ( init,  cond,  step,  body)

let cWhile ?(cond : target = [])
  ?(body : target = []) (_ : unit) : constr =
  let p_cond = cond in
  let p_body = body in
    Constr_while (p_cond, p_body)

let cDoWhile ?(body : target = [])
  ?(cond : target = [])  (_ : unit) : constr =
  let p_body = body in
  let p_cond = cond in
    Constr_do_while (p_cond, p_body)

let cIf ?(cond : target = [])
  ?(then_ : target = []) ?(else_ : target = []) (_ : unit) : constr =
  let p_cond = cond in
  let p_then = then_ in
  let p_else = else_ in
    Constr_if (p_cond, p_then, p_else)

let cThen : constr =
 Constr_target [cIf(); dThen]

let cOr (tgl : target list) : constr =
  Constr_or tgl

let cAnd (tgl : target list) : constr =
  Constr_and tgl

(* Converts a list of targets into a [target_list_pred] *)
let target_list_simpl (args : targets) : target_list_pred =
  let n = List.length args in
  make_target_list_pred
    (fun i -> if i < n then List.nth args i else [cStrict;bFalse]) (* TODO the else to name "target_none" *)
    (fun bs -> List.length bs = n && list_all_true bs)
    (fun () -> "target_list_simpl(" ^ (list_to_string (List.map target_to_string args) ^ ")"))

(* NOTE: the "_st" suffix means that the argument is a constraint and not a target
   --we might revisit this convention later if we find it not suitable *)

(* Converts a target into a [target_list_pred] that checks that at least one of the items in the list satisfies the given constraint *)
let target_list_one_st (tg : target) : target_list_pred =(* LATER: KEEP ONLY THIS *)
  make_target_list_pred
    (fun _i -> tg)
    (fun bs -> List.mem true bs)
    (fun () -> "target_list_one_st(" ^ (target_to_string tg) ^ ")")


(* Converts a target into a [target_list_pred] that checks that at least all the items in the list satisfies the given constraint *)
let target_list_all_st (tg : target) : target_list_pred = (* LATER: KEEP ONLY THIS *)
  make_target_list_pred
    (fun _i -> tg)
    (fun bs -> List.for_all (fun b -> b = true) bs)
    (fun () -> "target_list_all_st(" ^ (target_to_string tg) ^ ")")

(* Predicate that matches any list of arguments *)
let target_list_pred_default : target_list_pred =
  make_target_list_pred
    (fun _i -> [])
    list_all_true
    (fun () -> "target_list_pred_default")

(* [combine_args args args_pred] takes [args] as a [target_list_simpl] if it is nonempty,
   and else returns [args_pred]; raise an error if the two arguments have non-default values *)
let combine_args (args:targets) (args_pred:target_list_pred) : target_list_pred =
  match args with
  | [] -> args_pred
  | _ ->
      if args_pred != target_list_pred_default
        then fail None "cFunDef: can't provide both args and args_pred";
      target_list_simpl args


(* by default an empty name is no name *)
let cFunDef ?(args : targets = []) ?(args_pred : target_list_pred = target_list_pred_default) ?(body : target = []) ?(regexp : bool = false) (name : string) : constr =
  let ro = string_to_rexp_opt regexp false name TrmKind_Expr in
  Constr_decl_fun (ro, combine_args args args_pred, body)

(* toplevel fun declaration *)
let cTopFunDef
  ?(args : targets = []) ?(args_pred : target_list_pred = target_list_pred_default)
  ?(body : target = []) (name : string) : constr =
  cChain [ dRoot; cFunDef ~args ~args_pred ~body name ]

let cTypDef
  ?(substr : bool = false) ?(regexp : bool = false) (name : string) : constr =
  let ro = string_to_rexp_opt regexp substr name TrmKind_Expr in
  Constr_decl_type ro

let cDef (name : string) : constr =
  cOr [[cFunDef name];[cVarDef name];[cTypDef name]]

let cEnum ?(name : string = "")
  ?(substr : bool = false) ?(constants : (string * (target)) list = [])
  ?(regexp : bool = false) (_ : unit) : constr =
  let c_n = string_to_rexp_opt regexp substr name TrmKind_Expr in
  let cec_o =
    match constants with
    | [] -> None
    | _ ->
        let cec =
          List.map
            (fun (n, pl) -> (string_to_rexp_opt regexp substr n TrmKind_Expr, pl))
            constants
        in
        Some cec
  in
  Constr_decl_enum (c_n, cec_o)

let cSeq ?(args : targets = []) ?(args_pred:target_list_pred = target_list_pred_default) (_ : unit) : constr =
  Constr_seq (combine_args args args_pred)

let cVar ?(regexp : bool = false) ?(trmkind : trm_kind = TrmKind_Expr) (name : string) : constr =
  let ro = string_to_rexp_opt regexp false name trmkind in
  Constr_var ro

let cBool (b : bool) : constr =
    Constr_lit (Some (Lit_bool b))

let cInt (n : int) : constr =
    Constr_lit (Some (Lit_int n))

let cDouble (f : float) : constr =
    Constr_lit (Some (Lit_double f))

let cString (s : string) : constr =
    Constr_lit (Some (Lit_string s))

let cLit : constr =
   Constr_lit None

(* [cCall] can match all kind of function calls *)
let cCall ?(fun_  : target = []) ?(args : targets = []) ?(args_pred:target_list_pred = target_list_pred_default) ?(accept_encoded : bool = false) (name:string) : constr =
  let exception Argument_Error of string in
  let p_fun = match fun_ with
  | [] -> [cVar name]
  | _ ->
    begin match name with
    | "" -> fun_
    | _ -> raise (Argument_Error "Can't provide both the path and the name of the function")
    end in
  Constr_app (p_fun, combine_args args args_pred, accept_encoded)

(* [cFun] matches a function by its name; it cannot match primitive functions *)
let cFun ?(fun_  : target = []) ?(args : targets = []) ?(args_pred:target_list_pred = target_list_pred_default) (name:string) : constr =
  cCall ~fun_ ~args ~args_pred ~accept_encoded:false name

(* [cPrim] matches only primitive functions; use [cPrimFun] for matching primitive function calls. *)
let cPrim (p : prim) : constr =
  Constr_prim p

(* [cPrimFun ~args ~args_pred  p] matches only primitive function calls*)
let cPrimFun ?(args : targets = []) ?(args_pred:target_list_pred = target_list_pred_default) (p:prim) : constr =
   cCall ~fun_:[cPrim p] ~args ~args_pred ""

(* [cSet ~lhs ~rhs ()] matches set operations with left hand side [lhs] and right hand side [rhs], if right(left) hand side are
    left empty, then no contraint on the side of the set operation will be applied.
*)
let cSet ?(lhs : target = []) ?(rhs : target = []) (_ : unit) : constr =
  cPrimFun ~args:[lhs; rhs] (Prim_binop Binop_set)

let cGet ?(arg : target = []) () : constr = 
  cPrimFun ~args:[arg] (Prim_unop Unop_get)

(* [cMark m] matches all the ast nodes with annotation Mark m*)
let cMark (m : mark) : constr =
  Constr_mark ((fun m1 -> m1 = m), "exactly:" ^ m)

(* [cMarks ms] matches all the ast nodes with annotation Mark m when m is an element of ms *)
let cMarks (ms : mark list) : constr =
  Constr_mark ((fun m1 -> List.mem m1 ms), "one of:" ^ (Tools.list_to_string ms))

let cMarkSt (pred : mark -> bool) : constr =
  Constr_mark (pred, "such_that:" ^ "todo") 

(* [cMarkAny] matches all the ast nodes with annotation Mark m, where m can be any positive integer *)
let cMarkAny : constr =
  Constr_mark ((fun _ -> true), "any_mark")

(* [cLabel ~substr ~body ~regexp label] matches C labels*)
let cLabel ?(substr : bool = false) ?(body : target = []) ?(regexp : bool = false) (label : string) : constr =
  let ro = string_to_rexp_opt regexp substr label TrmKind_Expr in
  let p_body = body in
  Constr_label (ro, p_body)

let cGoto ?(label : string = "")
  ?(substr : bool = false) ?(regexp : bool = false) (_ : unit) : constr =
  let ro = string_to_rexp_opt regexp substr label TrmKind_Expr in
  Constr_goto ro

let cReturn_target ?(res : target = [])
  (_ : unit) : constr =
  let p_res =  res in
  Constr_return p_res

let cAbrtAny : abort_kind = Any

let cAbrtRet : abort_kind = Return

let cAbrtBrk : abort_kind = Break

let cAbrtCtn : abort_kind = Continue

let cAbort ?(kind : abort_kind = Any)
  (_ : unit) : constr =
  Constr_abort kind

let cReturn : constr =
  Constr_abort (cAbrtRet)

let cBreak : constr =
  Constr_abort (cAbrtBrk)

let cContinue : constr =
  Constr_abort (cAbrtCtn)
(*
  the empty list is interpreted as no constraint on the accesses
  accesses are reversed so that users give constraints on what they see
  *)
let cAccesses ?(base : target = [])
  ?(accesses : constr_access list = []) (_ : unit) : constr =
  let p_base =  base in
  let accesses =
    match accesses with | [] -> None | cal -> Some (List.rev cal)
  in
    Constr_access (p_base, accesses)

let cIndex ?(index : target = []) (_ : unit) : constr_access =
  let p_index =  index in
  Array_access p_index

let cField ?(field : string = "") ?(substr : bool = false) ?(regexp : bool = false)
  (_ : unit) : constr_access =
  let ro = string_to_rexp_opt regexp substr field TrmKind_Expr in
  Struct_access ro

let cAccess : constr_access = 
  Any_access

(* the empty list is interpreted as no constraint on the cases *)
let cSwitch ?(cond : target = [])
  ?(cases : (case_kind * (target)) list = []) (_ : unit) : constr =
  let p_cond =  cond in
  let c_cases =
    match cases with
    | [] -> None
    | _ -> Some (List.map (fun (k, pl) -> (k,  pl)) cases)
  in
    Constr_switch (p_cond, c_cases)

let cCase ?(value : target = []) (_ : unit) : case_kind =
  match value with
  | [] -> Case_any
  | _ -> Case_val ( value)

let cDefault : case_kind = Case_default

let dLHS : constr =
  cChain [cSet(); dArg 0]

let dRHS : constr =
  cChain [cSet (); dArg 1]

let cTargetInDepth (tg : target) : constr =
  Constr_target (Constr_depth DepthAny :: tg)





(******************************************************************************)
(*                          Target resolution                                 *)
(******************************************************************************)

(* NOW INCLUDED
let resolve_target = Constr.resolve_target
let resolve_target_between = Constr.resolve_target_between
*)


(******************************************************************************)
(*                          Apply on target operations                        *)
(******************************************************************************)

(* Type of transformations *)
module Transfo = struct
  type t = target -> unit
  type local = trm -> path -> trm
  type local_between = int -> local
end


let apply_on_path = Path.apply_on_path

(* [applyi_on_targets ~replace tr tg]: Apply a specific Generic over a target or a list of targets, keep track over the index of the target
      params:
        tg : target
        tr : transformation to be applied
      return:
        unit
*)

(* [applyi_on_transformed_targets transformer tr tg]: Apply a transformation [tr] on target [tg]
      params:
        transformer: change the resolved path so that more information about the context of the node is given
        tr: transformation to be applied at the nodes corresponedt to target [tg]
        tg: target
      return:
        unit
*)
let applyi_on_transformed_targets (transformer : path -> 'a) (tr : int -> trm -> 'a -> trm) (tg : target) : unit =
  Trace.apply (fun t ->
    let ps = resolve_target tg t in
    let marks = List.map (fun _ -> Mark.next()) ps in
    let t = List.fold_left2 (fun t p m -> apply_on_path (trm_add_mark m) t p) t ps marks in
    Tools.foldi( fun imark t m ->
      match resolve_target [cMark m] t with
      | [] -> fail None "applyi_on_transformed_targets: a mark disappeared"
      | [p] -> let t = apply_on_path (trm_remove_mark m) t p in
        tr imark t (transformer p)
      | _ -> fail None "applyi_on_transformed_targets: a mark was duplicated"
    ) t marks)

(* [apply_on_transformed_targets ~replace_top transformer tr tg]:
   Same as [apply_to_transformed_targets] except that there is some processing performed on each of the explicit path.
    params:
      transformer: ..
      tr: transformation to be applied
      tg: target
    return:
      unit
*)

let apply_on_transformed_targets (transformer : path -> 'a) (tr : 'a -> trm -> trm) (tg : target) : unit =
  applyi_on_transformed_targets  transformer (fun _i t descr -> tr descr t) tg

(* [applyi_on_transformed_targets tr tg]: similar as applyi_on_transformed_targets except that this function does
    note modfiy the resolved path
      params:
        tr: transformation
        tg: target
*)
let applyi_on_targets (tr : int -> trm -> path -> trm) (tg : target) : unit =
  applyi_on_transformed_targets (fun p -> p) tr tg




(* [apply_on_targets ~replace tr tg]: esentially the same as applyi_on_targets, but without keeping track over the index of the target
      params:
        tg : target
        tr : transformation to be applied
      return:
        unit
*)
let apply_on_targets (tr : trm -> path -> trm) (tg : target) : unit =
  applyi_on_targets (fun _i t dl -> tr t dl) tg

let applyi_on_transformed_targets_between (transformer : path * int -> 'a) (tr : int -> trm -> 'a -> trm) (tg : target) : unit =
  Trace.apply( fun t ->
  let ps = resolve_target_between tg t in
  let marks = List.map (fun _ -> Mark.next ()) ps in
  let t = List.fold_left2 (fun t (p_to_seq, i) m -> apply_on_path (trm_add_mark_between i m) t p_to_seq ) t ps marks in
  Tools.foldi (fun imark t m ->
    match resolve_target [cMark m] t with
    | [] -> fail None "applyi_on_transformed_targets_between: a mark disappeared"
    | [p_to_seq] ->
      let t_seq, _ = resolve_path p_to_seq t in
      let i = begin match get_mark_index m t_seq with | Some i -> i | None -> fail t_seq.loc "applyi_on_transformed_targets_between: could not get the between index" end in
      let t = apply_on_path (trm_remove_mark_between m) t p_to_seq in
      tr imark t (transformer (p_to_seq,i))
    | _ -> fail None "applyi_on_transformed_targets_between: a mark was duplicated"
  ) t marks)

let apply_on_transformed_targets_between (transformer: path * int -> 'a) (tr : trm -> 'a -> trm) (tg : target) : unit =
  applyi_on_transformed_targets_between transformer (fun _i t descr -> tr t descr) tg


(* [apply_on_targets_between ~replace_top tr tg]: Similar to apply_on_targets, but the function considers the index too
      params:
        tr : transformation to be applied
        tg : target
      return:
        unit
*)
let applyi_on_targets_between (tr : int -> trm -> path * int -> trm) (tg : target) : unit =
  applyi_on_transformed_targets_between (fun (p,i) -> (p,i)) tr tg

let apply_on_targets_between (tr : trm -> 'a -> trm) (tg : target) : unit =
  applyi_on_targets_between (fun _i t pk -> tr t pk) tg



(******************************************************************************)
(*                                   Show                                     *)
(******************************************************************************)

(* [target_show_aux id t]: adds an annotation [trm_decoration]
   carrying the information [id] around the term t.
*)
let target_show_aux (id : int) (t : trm) : trm =
  let show_mark = "show_mark " ^ (string_of_int id) in
  trm_add_mark show_mark t

(* [target_show_transfo id t p]: adds an annotation [trm_decoration]
   carrying the information [id] around the term at path [p] in the term [t]. *)
let target_show_transfo (id : int) : Transfo.local =
  apply_on_path (target_show_aux id)

(* [target_between_show_aux id k t]: adds a decorated semi-column with identifier [id]
   at position [k] in the sequence described by the term [t]. *)
let target_between_show_aux (id : int) (k : int) (t : trm) : trm =
    trm_add_mark_between k (string_of_int id) t 

(* [target_between_show_transfo id k t p]: adds a decorated semi-column with identifier [id]
   at position [k] in the sequence at path [p] in the term [t]. *)
let target_between_show_transfo (id : int) : Transfo.local_between =
  fun (k:int) -> apply_on_path (target_between_show_aux id k)

(* [show ~line:int tg] is a transformation for visualizing targets.
   The operation only executes if the command line argument [-exit-line]
   matches the [line] argument provided to the function. Otherwise, it is a noop.
   There is no need for a prefix [!!] or [!!!] to the front of the [show]
   function, because it is recognized as a special function by the preprocessor
   that generates the [foo_with_lines.ml] instrumented source. *)
let show ?(line : int = -1) ?(reparse : bool = true) (tg : target) : unit =
        if reparse then reparse_alias();
        let should_exit = (Flags.get_exit_line() = Some line) in
        if should_exit then begin 
          if Constr.is_target_between tg then begin
            applyi_on_targets_between (fun i t (p,k) -> 
              target_between_show_transfo i k t p) tg
              end
          else begin
            applyi_on_targets (fun i t p -> target_show_transfo i t p) tg
          end;
          dump_diff_and_exit()
        end

(* DEPRECATED *)
(* let show ?(line : int = -1) ?(reparse : bool = false) (tg : target) : unit =
  only_interactive_step line ~reparse (fun () ->
    if Constr.is_target_between tg then begin
      applyi_on_targets_between (fun i  t (p,k) ->
        target_between_show_transfo i k t p) tg
    end else begin
      applyi_on_targets (fun i t p ->
        target_show_transfo i t p) tg
    end) *)

(** [force_reparse_after tr] is a wrapper to invoke for forcing the reparsing
    after a transformation. For example because it modifies type definitions.
    See example in [Struct.inline]. *)
let force_reparse_after (tr : Transfo.t) : Transfo.t =
  fun (tg : target) ->
    tr tg;
    Trace.reparse()

