open Ast
open Ast_to_c
open Str
open Tools
(******************************************************************************)
(*                                  Path AST                                  *)
(******************************************************************************)


     
(* explicit path in trm ast = list of directions *)
(* todo: find better name (trm_path?) *)
type path = dir list

and dir =
  (* nth: go to nth element in seq, array, struct *)
  | Dir_nth of int
  (* cond: used for if, loops and switch *)
  | Dir_cond
  (* if *)
  | Dir_then
  | Dir_else
  (*
    body: used for loops, definitions, return, labelled terms or switch case
      -> directions for while loop: cond and body
   *)
  | Dir_body
  (* for *)
  | Dir_for_init
  | Dir_for_step
  (* app *)
  | Dir_app_fun
  (* arg for fun application and declaration *)
  | Dir_arg of int
  (* name of declared var/fun or label *)
  | Dir_name
  (*
    case group in switch
    Dir_case (n, d) = follow d in nth case group
   *)
  | Dir_case of int * case_dir
  (* constant in enum declaration *)
  | Dir_enum_const of int * enum_const_dir

and case_dir =
  | Case_name of int
  | Case_body

and enum_const_dir =
  | Enum_const_name
  | Enum_const_val

let dir_to_string (d : dir) : string =
  match d with
  | Dir_nth n -> "Dir_nth " ^ (string_of_int n)
  | Dir_cond -> "Dir_cond"
  | Dir_then -> "Dir_then"
  | Dir_else -> "Dir_else"
  | Dir_body -> "Dir_body"
  | Dir_for_init -> "Dir_for_init"
  | Dir_for_step -> "Dir_for_step"
  | Dir_app_fun -> "Dir_app_fun"
  | Dir_arg n -> "Dir_arg " ^ (string_of_int n)
  | Dir_name -> "Dir_name"
  | Dir_case (n, cd) ->
     let s_cd =
       match cd with
       | Case_name n -> "Case_name " ^ (string_of_int n)
       | Case_body -> "Case_body"
     in
     "Dir_case (" ^ (string_of_int n) ^ ", " ^ s_cd ^ ")"
  | Dir_enum_const (n, ecd) ->
     let s_ecd =
       match ecd with
       | Enum_const_name -> "Enum_const_name"
       | Enum_const_val -> "Enum_const_val"
     in
     "Dir_enum_const (" ^ (string_of_int n) ^ ", " ^ s_ecd ^ ")"


let list_to_string ?(sep:string=";") ?(bounds:string list = ["[";"]"])(l : string list) : string =
  let rec aux = function
    | [] -> ""
    | [s] -> s
    | s1 :: s2 :: sl -> s1 ^ sep ^ " " ^ aux (s2 :: sl)
  in
  (List.nth bounds 0) ^ aux l ^ (List.nth bounds 1)


let string_of_explicit_path (dl : path) : string =
  list_to_string (List.map dir_to_string dl)

(*
  comparison functions for path sorting
  the order between directions does not matter
  only constraint: when one path is the prefix of the other, it must be
  considered "greater"
 *)
let compare_dir (d : dir) (d' : dir) : int =
  match d, d' with
  | Dir_nth n, Dir_nth m -> compare n m
  | Dir_arg n, Dir_arg m -> compare n m
  | Dir_case (n, cd), Dir_case (m, cd') ->
     let cn = compare n m in
     if cn <> 0 then cn else
       begin match cd, cd' with
       | Case_name i, Case_name j -> compare i j
       | Case_body, Case_body -> 0
       | Case_name _, _ -> -1
       | _, Case_name _ -> 1
       end
  | Dir_enum_const (n, ecd), Dir_enum_const (m, ecd') ->
     let cn = compare n m in
     if cn <> 0 then cn else
       begin match ecd, ecd' with
       | d, d' when d = d' -> 0
       | Enum_const_name, _ -> -1
       | Enum_const_val, _ -> 1
       end
  | d, d' when d = d' -> 0
  | Dir_nth _, _ -> -1
  | _, Dir_nth _ -> 1
  | Dir_cond, _ -> -1
  | _, Dir_cond -> 1
  | Dir_then, _ -> -1
  | _, Dir_then -> 1
  | Dir_else, _ -> -1
  | _, Dir_else -> 1
  | Dir_body, _ -> -1
  | _, Dir_body -> 1
  | Dir_for_init, _ -> -1
  | _, Dir_for_init -> 1
  | Dir_for_step, _ -> -1
  | _, Dir_for_step -> 1
  | Dir_app_fun, _ -> -1
  | _, Dir_app_fun -> 1
  | Dir_arg _, _ -> -1
  | _, Dir_arg _ -> 1
  | Dir_name, _ -> -1
  | _, Dir_name -> 1
  | Dir_case _, _ -> -1
  | _, Dir_case _ -> 1

let rec compare_path (dl : path) (dl' : path) : int =
  match dl, dl' with
  | [], [] -> 0
  | [], _ -> 1
  | _, [] -> -1
  | d :: dl, d' :: dl' ->
     let cd = compare_dir d d' in
     if cd = 0 then compare_path dl dl' else cd

(*
  regexps with their string description, a boolean for exact/subexpression
  matching, and the kind of expression to match (any, instruction, or
  subexpression)
 *)
(* TODO: replace exact with standard name *)
type rexp = {desc : string; exp : regexp; exact : bool; only_instr : bool}

type paths = path list

(* Relative target to term *)
type target_relative =
    | TargetAt (* default value *)
    | TargetFirst
    | TargetLast
    | TargetBefore
    | TargetAfter

(* Number of targets to match *)
type target_occurences =
    | ExpectedOne  (* = 1 occurence, the default value *)
    | ExpectedNb of int (* exactly n occurrences *)
    | ExpectedMulti    (* > 0 number of occurrences *)
    | ExpectedAnyNb  (* any number of occurrences *)

(* a target is a list of constraints to identify subterms *)
type target = constr list

and constr =
  (*
    target constraints:
    - directions to follow
    - list constraint: a constraint is matched against all elements of the list
      used for seq elements and function arguments
      the user provides a function that, given the results of this constraint
      matching, returns the ranks of the elements to explore next
    - include file to explore
   *)
  | Constr_dir of dir
  | Constr_list of target * (bool list -> int list)
  | Constr_include of string
  (*
    matching constraint: match against regexp
    the user may also match against a string through smart constructors
   *)
  (* todo: beware of multiline matching *)
  | Constr_regexp of rexp
  (*
    node related constraints (constraints are expressed using a target):
    - for loop: constraints on the initialisation, the condition, the step
    instruction and the body
    - while loop: constraints on the condition and the body
    - if statement: constraints on the condition and each branch
    - decl_var: constraints on the name and on the initialisation
    - decl_fun: constraints on the name, the arguments and the body
    - decl_type: constraint on the name
    - decl_enum: constraints on the name and the constants
    - seq: constraint matched against each all elements of a seq
      depending on the result, a function decides whether the constraint holds
    - var: constraint on the name
    - lit: literal
    - app: constraints on the function and on the list of arguments
    - label: constraints on the label and on the term
    - goto: constraint on the label
    - return: constraint on the returned value
    - abort: matches abort nodes (return, break and continue)
    - accesses: constraints on a succession of accesses + on the base
    - switch: constraints on the condition and on the cases
   *)
  (* for: init, cond, step, body *)
  | Constr_for of target * target * target * target
  (* while: cond, body *)
  | Constr_while of target * target
  (* if: cond, then, else *)
  | Constr_if of target * target * target
  (* decl_var: name, body *)
  | Constr_decl_var of constr_name * target
  (* decl_fun: name, args, body *)
  | Constr_decl_fun of constr_name * constr_list * target
  (* decl_type: name *)
  | Constr_decl_type of constr_name
  (* decl_enum: name, constants *)
  | Constr_decl_enum of constr_name * constr_enum_const
  (* seq *)
  | Constr_seq of constr_list
  (* var *)
  | Constr_var of constr_name
  (* lit *)
  | Constr_lit of lit
  (* app: function, arguments *)
  | Constr_app of target * constr_list
  (* label *)
  | Constr_label of constr_name * target
  (* goto *)
  | Constr_goto of constr_name
  (* return *)
  | Constr_return of target
  (* abort *)
  | Constr_abort of abort_kind
  (* accesses: base, accesses *)
  | Constr_access of target * constr_accesses
  (* switch: cond, cases *)
  | Constr_switch of target * constr_cases
  (* TODO: Constraint for types? *)
  (* Target relative to another trm *)
  | ConstrRelative of target_relative
  (* Number of  occurrences expected  *)
  | ConstrOccurences of target_occurences
  (* List of constraints *)
  | ConstrChain of constr list

and constr_name = rexp option

and constr_list = target * (bool list -> bool)

and constr_enum_const = ((constr_name * target) list) option

and constr_accesses = (constr_access list) option

and constr_access =
  (* array indices may be arbitrary terms *)
  | Array_access of target
  (* struct fields are strings *)
  | Struct_access of constr_name
  | Any_access

(* for each case, its kind and a constraint on its body *)
and constr_cases = ((case_kind * target) list) option

and case_kind =
  (* case: value *)
  | Case_val of target
  | Case_default
  | Case_any

and abort_kind =
  | Any
  | Return
  | Break
  | Continue


type target_simple = target (* Without ConstrRelative, ConstrOccurences, ConstrChain *)

(* Advanced path search *)
type target_struct = {
   target_path: target_simple; (* this path contains no cMulti/cNb/cBefore/etc.., only cStrict can be there *)
   target_relative: target_relative;
   target_occurences: target_occurences; }

(* Flatten all the constrainst of type ConstrChain *)
let target_flatten(tg : target) : target =
    let rec aux cl = match cl with 
    | [] -> []
    | x :: xs -> begin match x with 
                | ConstrChain tr ->  tr @ aux xs
                | _ -> x :: aux xs 
                end 
    in aux tg
                     
(* Convert a target into a target struct  *)
let target_to_target_struct(tg : target) : target_struct =
  let tg = target_flatten tg in 
  let relative = ref None in 
  let occurences = ref None in
  let process_constr (c : constr) : unit =
    match c with 
    | ConstrRelative tr ->
      begin match !relative with 
      | None -> relative := Some tr;
      | Some _ -> fail None  "ConstrRelative provided twice in path"
      end 
    | ConstrOccurences oc -> 
      begin match !occurences with
      | None -> occurences := Some oc;
      | _ -> fail None "ConstrOccurences provide twice in path"
      end
    | _ -> ()
   (* Check if relative constraint are applied once and the number of occurences is unique *)
   in List.iter process_constr tg;
   
   (* Return a target_struct *)
   {   target_path = List.filter (function | ConstrRelative _ | ConstrOccurences _ -> false | _ -> true) tg;
       target_relative = begin match !relative with | None -> TargetAt | Some tg -> tg end;
       target_occurences = begin match !occurences with | None -> ExpectedOne | Some oc -> oc end; 
   }

let regexp_to_string (r : rexp) : string =
  (if r.exact then "Exact " else "Sub ") ^
    (if r.only_instr then "OnlyInstr \"" else "\"") ^ r.desc ^ "\""

let rec constraint_to_string (c : constr) : string =
  match c with
  | Constr_dir d -> dir_to_string d
  | Constr_list (p_elt, _) -> "List (" ^ target_to_string p_elt ^ ")"
  | Constr_include s -> "Include " ^ s
  | Constr_regexp r -> "Regexp " ^ regexp_to_string r
  | Constr_for (p_init, p_cond, p_step, p_body) ->
     let s_init = target_to_string p_init in
     let s_cond = target_to_string p_cond in
     let s_step = target_to_string p_step in
     let s_body = target_to_string p_body in
     "For (" ^ s_init ^ ", " ^ s_cond ^ ", " ^ s_step ^ ", " ^ s_body ^ ")"
  | Constr_while (p_cond, p_body) ->
     let s_cond = target_to_string p_cond in
     let s_body = target_to_string p_body in
     "While (" ^ s_cond ^ ", " ^ s_body ^ ")"
  | Constr_if (p_cond, p_then, p_else) ->
     let s_cond = target_to_string p_cond in
     let s_then = target_to_string p_then in
     let s_else = target_to_string p_else in
     "If (" ^ s_cond ^ ", " ^ s_then ^ ", " ^ s_else ^ ")"
  | Constr_decl_var (name, p_body) ->
     let s_name =
       match name with | None -> "_" | Some r -> regexp_to_string r
     in
     let s_body = target_to_string p_body in
     "Decl_var (" ^ s_name ^ ", " ^ s_body ^ ")"
  | Constr_decl_fun (name, (p_args, _), p_body) ->
     let s_name =
       match name with | None -> "_" | Some r -> regexp_to_string r
     in
     let s_args = target_to_string p_args in
     let s_body = target_to_string p_body in
     "Decl_fun (" ^ s_name ^ ", " ^ s_args ^ ", " ^ s_body ^ ")"
  | Constr_decl_type name ->
     let s_name =
       match name with | None -> "_" | Some r -> regexp_to_string r
     in
     "Decl_type " ^ s_name
  | Constr_decl_enum (name, c_const) ->
     let s_name =
       match name with | None -> "_" | Some r -> regexp_to_string r
     in
     let s_const =
       match c_const with
       | None -> "_"
       | Some np_l ->
          let sl =
            List.map
              (fun (n, p) ->
                let s_n =
                  match n with | None -> "_" | Some r -> regexp_to_string r
                in
                let s_p = target_to_string p in
                "(" ^ s_n ^ ", " ^ s_p ^ ")"
              )
              np_l
          in
          list_to_string sl
     in
     "Decl_enum (" ^ s_name ^ ", " ^ s_const ^ ")"
  | Constr_seq (p_elt, _) ->
     let s = target_to_string p_elt in
     "Seq (" ^ s ^ ")"
  | Constr_var name ->
     "Var " ^ (match name with | None -> "_" | Some r -> regexp_to_string r)
  | Constr_lit l ->
     let s =
       begin match l with
       | Lit_unit ->  "()"
       | Lit_uninitialized -> "?"
       | Lit_bool b -> string_of_bool b
       | Lit_int n -> string_of_int n
       | Lit_double d -> string_of_float d
       | Lit_string s -> s
       end
     in "Lit " ^ s
  | Constr_app (p_fun, (p_args, _)) ->
     let s_fun = target_to_string p_fun in
     let s_args = target_to_string p_args in
     "App (" ^ s_fun ^ ", " ^ s_args ^ ")"
  | Constr_label (so, p_body) ->
     let s_label =
       match so with | None -> "_" | Some r -> regexp_to_string r
     in
     let s_body = target_to_string p_body in
     "Label (" ^ s_label ^ ", " ^ s_body ^ ")"
  | Constr_goto so ->
     let s_label =
       match so with | None -> "_" | Some r -> regexp_to_string r
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
  | Constr_access (p_base, ca) ->
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
  | ConstrRelative tr ->
    begin match tr with 
    | TargetAt -> "TargetAt"
    | TargetFirst -> "TargetFirst"
    | TargetLast -> "TargetLast"
    | TargetBefore -> "TargetBefore"
    | TargetAfter -> "TargetAfter"
    end 
  | ConstrOccurences oc ->
    begin match oc with 
    | ExpectedOne -> "ExpectedOne"
    | ExpectedNb n-> "ExpectedNb " ^ string_of_int(n)
    | ExpectedMulti -> "ExpectedMulti"
    | ExpectedAnyNb -> "ExpectedAnyNb"
    end
  | ConstrChain cl ->
    let string_cl = List.map constraint_to_string cl in 
    list_to_string string_cl


and target_to_string (tg : target) : string =
  list_to_string (List.map constraint_to_string tg)

and access_to_string (ca : constr_access) : string =
  match ca with
  | Array_access p_index ->
     let s_index = target_to_string p_index in
     "Array_access " ^ s_index
  | Struct_access so ->
     let s_field =
       match so with | None -> "_" | Some r -> regexp_to_string r
     in
     "Struct_access " ^ s_field
  | Any_access -> "Any_access"

(******************************************************************************)
(*                        Smart constructors for targets                        *)
(******************************************************************************)

(* TODO: Remove all the occurrences of List.flatten they are not needed anymore*)
let (@@) (l:'a list) = List.flatten l     

module Path_constructors =
  struct
    (*
      a smart constructor builds a target
      thus, the user provides a target list using them
      this list is then flattened to call resolve_target
      unit args are used because of optional arguments
     *)
    let cBefore : constr = 
      ConstrRelative TargetBefore

    let cAfter : constr =
      ConstrRelative TargetAfter

    let cFirst : constr =
      ConstrRelative TargetFirst

    let cLast : constr = 
      ConstrRelative TargetLast

    let cMulti : constr =
      ConstrOccurences ExpectedMulti
    
    let cAnyNb : constr =
       ConstrOccurences ExpectedAnyNb
    
    let cNb (nb : int) : constr =
       ConstrOccurences (ExpectedNb nb)
    
    (* directions *)
    let cNth (n : int) : target =
       [Constr_dir (Dir_nth n)]
    let cCond (_ : unit) : target =
       [Constr_dir Dir_cond]
    let cThen (_ : unit) : target =
       [Constr_dir Dir_then]
    let cElse (_ : unit) : target =
       [Constr_dir Dir_else]
    let cBody (_ : unit) : target =
       [Constr_dir Dir_body]
    let cInit (_ : unit) : target =
       [Constr_dir Dir_for_init]
    let cStep (_ : unit) : target =
       [Constr_dir Dir_for_step]
    let cAppFun (_ : unit) : target =
       [Constr_dir Dir_app_fun]
    let cArg (n : int) : target =
       [Constr_dir (Dir_arg n)]
    let cName (_ : unit) : target =
       [Constr_dir Dir_name]
    let cDirCase (n : int) (cd : case_dir) : target =
       [Constr_dir (Dir_case (n, cd))]
    let cCaseName (n : int) : case_dir = Case_name n
    let cCaseBody : case_dir = Case_body
    let cEnumConst (n : int)
      (ecd : enum_const_dir) : target =
       [Constr_dir (Dir_enum_const (n, ecd))]
    let cEnumConstName : enum_const_dir = Enum_const_name
    let cEnumConstVal : enum_const_dir = Enum_const_val

    (* constraints *)
    let cList_int (tgl : target list)
      (next : bool list -> int list) : target =
       [Constr_list ((@@) tgl, next)]

    (* allow to use boolean functions *)
    let cList (tgl : target list)
      (next : bool list -> bool list) : target =
      let rec int_of_bool_list (n : int) = function
        | [] -> []
        | b :: bl ->
           let il = int_of_bool_list (n + 1) bl in
           if b then n :: il else il
      in
      
        [Constr_list ((@@) tgl, fun bl -> int_of_bool_list 0 (next bl))]

    (* continue on the first matching instruction *)
    let cFirst (tgl : target list) : target =
      let rec next = function
        | [] -> []
        | true :: bl -> true :: List.map (fun _ -> false) bl
        | false :: bl -> false :: next bl
      in
      cList tgl next

    (* after operator *)
    (*
      after
     *)
    let rec after_bool (bl : bool list) : bool list =
      match bl with
      | [] -> []
      | [_] -> [false]
      | false :: bl -> false :: after_bool  bl
      | true :: _ :: bl ->
         let bl' = List.map (fun _ -> true) bl
         in
         false :: true :: bl'


    let before_aux (bl : bool list) : int list =
      match get_index true bl with
      | None -> []
      | Some 0 -> []
      | Some n -> List.init n (fun m -> m)

    let cInclude (s : string) : target =
       [Constr_include s]

    let rexp_of_string ?(only_instr : bool = true) ?(exact : bool = true)
          (s : string) : rexp =
      let exp = if exact then s ^ "$" else s in
      {desc = s; exp = Str.regexp exp; exact; only_instr}

    let cRegexp ?(exact : bool = true)
      ?(only_instr : bool = true) (s : string) : target =
       [Constr_regexp (rexp_of_string ~only_instr ~exact s)]

    (* exactly match the string/regexp described by s *)
    let cStr ?(regexp : bool = false)
      (s : string) : target =
      cRegexp ~only_instr:false (if regexp then s else Str.quote s)

    (*
      match the string/regexp only on instructions
      by default this is not an exact match
     *)
    let cInstrSubstr ?(exact : bool = false)
      ?(regexp : bool = false) (s : string) : target =
      cRegexp  ~exact ~only_instr:true
        (if regexp then s else Str.quote s)
    let rexp_opt_of_string ?(exact : bool = true) (s : string) : rexp option =
      if s = "" then None else Some (rexp_of_string ~only_instr:false ~exact s)

    let cVarDef ?(name : string = "")
      ?(exact : bool = true) ?(body : target list = []) (_ : unit) : target =
      let ro = rexp_opt_of_string ~exact name in
      let p_body = (@@) body in
       [Constr_decl_var (ro, p_body)]

    let cFor ?(init : target list = [])
      ?(cond : target list = []) ?(step : target list = []) ?(body : target list = []) ?(name : string = "")
      (_ : unit) : target =
      let init =
         match name, init with
         | "",[] -> init (*failwith "cFor: Need to provide the name or init"*)
         | "", _ -> init
         | _, [] -> [cVarDef ~name ()]
         | _, _::_ -> failwith "cFor: cannot provide both name and init"
         in
      
       [Constr_for ((@@) init, (@@) cond, (@@) step, (@@) body)]

    let cWhile ?(cond : target list = [])
      ?(body : target list = []) (_ : unit) : target =
      let p_cond = (@@) cond in
      let p_body = (@@) body in
       [Constr_while (p_cond, p_body)]

    let cIf ?(cond : target list = [])
      ?(then_ : target list = []) ?(else_ : target list = []) (_ : unit) : target =
      let p_cond = (@@) cond in
      let p_then = (@@) then_ in
      let p_else = (@@) else_ in
       [Constr_if (p_cond, p_then, p_else)]

    (* by default an empty name is no name *)
    

    let cFun ?(name : string = "")
      ?(exact : bool = true) ?(args : target list = [])
      ?(validate : bool list -> bool = fun _ -> true) ?(body : target list = [])
      (_ : unit) : target =
      let ro = rexp_opt_of_string ~exact name in
      let p_args = (@@) args in
      let p_body = (@@) body in
       [Constr_decl_fun (ro, (p_args, validate), p_body)]

    (* toplevel fun declaration *)
    let cTopFun ?(name : string = "") ?(exact : bool = true)
      ?(args : target list = []) ?(validate : bool list -> bool = fun _ -> true)
      ?(body : target list = []) (_ : unit) : target =
      (*
        structure of toplevel term:
        seq (del_instr)
          [
            seq (no_brace)
              [
                ... include ...
                t_top
              ]
            ... del_instr ...
          ]
        thus:
          1. find the list that contains the list that contains the function and
             explore the inner list
          2. find the function
       *)
      (cList 
         [cList 
           [cFun  ~name ~exact ~args ~validate ~body ()]
           (fun l -> l)
         ]
         (fun l -> l)
      ) ++
      (cList 
         [cFun  ~name ~exact ~args ~validate ~body ()]
         (fun l -> l)
      )

    let cType ?(name : string = "")
      ?(exact : bool = true) (_ : unit) : target =
      let ro = rexp_opt_of_string ~exact name in
       [Constr_decl_type ro]

    let cEnum ?(name : string = "")
      ?(exact : bool = true) ?(constants : (string * (target list)) list = [])
      (_ : unit) : target =
      let c_n = rexp_opt_of_string ~exact name in
      let cec_o =
        match constants with
        | [] -> None
        | _ ->
           let cec =
             List.map
               (fun (n, pl) -> (rexp_opt_of_string ~exact n, (@@) pl))
               constants
           in
           Some cec
      in
       [Constr_decl_enum (c_n, cec_o)]

    let cSeq ?(args : target list = [])
      ?(validate : bool list -> bool = fun _ -> true) (_ : unit) : target =
      let p_args = (@@) args in
       [Constr_seq (p_args, validate)]

    let cVar ?(name : string = "")
      ?(exact : bool = true) (_ : unit) : target =
      let ro = rexp_opt_of_string ~exact name in
       [Constr_var ro]

    let cBool (b : bool) : target =
       [Constr_lit (Lit_bool b)]
    let cInt (n : int) : target =
       [Constr_lit (Lit_int n)]
    let cDouble (f : float) : target =
       [Constr_lit (Lit_double f)]
    let cString (s : string) : target =
       [Constr_lit (Lit_string s)]
    let cPrim (p : prim) : target =
      cStr  (ast_to_string (trm_prim p))

    let cApp ?(name : string = "")
      ?(fun_ : target list = []) ?(args : target list = [])
      ?(validate : bool list -> bool = fun _ -> true) (_ : unit) : target =
      let exception Argument_Error  of string in  
      let p_fun =
      match name, fun_ with 
      | "",_ -> (@@) fun_ 
      | _, [] -> cVar  ~name ()
      | _,_ ->  raise (Argument_Error "Can't provide both the target and the name of the function")
      
      
      in
      let p_args = (@@) args in
       [Constr_app (p_fun, (p_args, validate))]

    let cLabel ?(label : string = "")
      ?(exact : bool = true) ?(body : target list = []) (_ : unit) : target =
      let ro = rexp_opt_of_string ~exact label in
      let p_body = (@@) body in
       [Constr_label (ro, p_body)]

    let cGoto ?(label : string = "")
      ?(exact : bool = true) (_ : unit) : target =
      let ro = rexp_opt_of_string ~exact label in
       [Constr_goto ro]

    let cReturn ?(res : target list = [])
      (_ : unit) : target =
      let p_res = (@@) res in
       [Constr_return p_res]

    let cAbort ?(kind : abort_kind = Any)
      (_ : unit) : target =
       [Constr_abort kind]

    let cAbrtAny : abort_kind = Any
    let cAbrtRet : abort_kind = Return
    let cAbrtBrk : abort_kind = Break
    let cAbrtCtn : abort_kind = Continue

    (*
      the empty list is interpreted as no constraint on the accesses
      accesses are reversed so that users give constraints on what they see
     *)
    let cAccesses ?(base : target list = [])
      ?(accesses : constr_access list = []) (_ : unit) : target =
      let p_base = (@@) base in
      let accesses =
        match accesses with | [] -> None | cal -> Some (List.rev cal)
      in
       [Constr_access (p_base, accesses)]

    let cIndex ?(index : target list = []) (_ : unit) : constr_access =
      let p_index = (@@) index in
      Array_access p_index

    let cField ?(field : string = "") ?(exact : bool = true)
      (_ : unit) : constr_access =
      let ro = rexp_opt_of_string ~exact field in
      Struct_access ro

    let cAccess : constr_access = Any_access

    (* the empty list is interpreted as no constraint on the cases *)
    let cSwitch ?(cond : target list = [])
      ?(cases : (case_kind * (target list)) list = []) (_ : unit) : target =
      let p_cond = (@@) cond in
      let c_cases =
        match cases with
        | [] -> None
        | _ -> Some (List.map (fun (k, pl) -> (k, (@@) pl)) cases)
      in
       [Constr_switch (p_cond, c_cases)]

    let cCase ?(value : target list = []) (_ : unit) : case_kind =
      match value with
      | [] -> Case_any
      | _ -> Case_val ((@@) value)

    let cDefault : case_kind = Case_default

    let cSet ?(lhs : target list = [])
      ?(rhs : target list = []) (_ : unit) : target =
      (@@)
        [
          (* first check that lhs is the first of two arguments *)
          cApp  ~args:lhs
            ~validate:(function | true :: [_] -> true | _ -> false) ();
          cApp  ~args:rhs
            ~validate:(function | _ :: [true] -> true | _ -> false) ();
          (* finally check that the function prints as "=" *)
          cApp  ~fun_:[cStr  "="] ()
        ]
  end

(******************************************************************************)
(*                              Target resolution                               *)
(******************************************************************************)

(*
  Particular case for target resolution: heap allocated variables
  Patterns:
    - declaration: seq annotated with Heap_allocated containing decl +
      optional initialisation annotated with Initialisation_instruction
    - usage: particular use of get/access because variables are pointers
      in practice, only dereferencing has to be taken into account: array and
      struct get/access are better matched with regexp
    - elimination:
      + in return instruction: seq annotated with Delete_instructions containing
        delete instructions annotated with Heap_allocated + abort instruction
      + at the end of scopes: seq annotated with Delete_instructions containing
        last instruction of the scope + annotated delete instructions
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
let forget_heap_alloc (t : trm) : trm =
  let loc = t.loc in
  match t.annot with
  | Some Delete_instructions ->
     (*
          t is either a sequence of delete instructions followed by a return
          or an instruction followed by a sequence of delete instructions
          in both case we do as if there was only the single instruction
      *)
     begin match t.desc with
     | Trm_seq tl ->
        (* all delete instructions are annotated by Heap_allocated *)
        begin match filter_out_heap_alloc tl with
        | [t'] -> t'
        | _ -> fail loc "forget_heap_alloc: bad delete list"
        end
     | _ -> fail loc "forget_heap_alloc: delete instructions should form a list"
     end
  | Some Heap_allocated ->
     (*
          t is either a sequence heap allocating a variable
          or the dereferencing of such a variable
          in both cases we replace the pointer with the value it points to
      *)
     begin match t.desc with
     (* first case: dereferencing, the user only sees the variable (t') *)
     | Trm_apps (_, [t']) -> t'
     (* second case: declaration *)
     | Trm_seq tl ->
        (* two subcases: with or without initial value *)
        begin match tl with
        (* first subcase: no initial value (init = new …) *)
        | [{desc = Trm_decl (Def_var ((x, {ty_desc = Typ_ptr tx; _}), _));
            _}] ->
           trm_decl ~loc ~is_instr:true
             (Def_var ((x, tx), trm_lit ~loc Lit_uninitialized))
        (* second subcase: initialisation *)
        | [{desc = Trm_decl (Def_var ((x, {ty_desc = Typ_ptr tx; _}), _)); _};
           {desc = Trm_apps (_, [{desc = Trm_var y; _}; init]); _}]
             when y = x ->
           trm_decl ~loc ~is_instr:true (Def_var ((x, tx), init))
        | _ -> fail loc "forget_heap_alloc: bad heap allocation"
        end
     | _ -> fail loc "forget_heap_alloc: bad heap_alloc instruction"
     end
  | _ -> fail loc "forget_heap_alloc: no heap_alloc in term"

(* return the last element of a list together with its index *)
let last (l : 'a list) : int * 'a =
  let rec aux n = function
    | [] -> failwith "last: empty list"
    | [a] -> (n, a)
    | _ :: b :: al -> aux (n + 1) (b :: al)
  in
  aux 0 l

(* applies a continuation to the nth element of l if it exists *)
let app_to_nth (loc : location) (l : 'a list) (n : int) (cont : 'a -> 'b) : 'b =
  try
    match List.nth_opt l n with
    | None ->
       fail loc
         ("app_to_nth: not enough elements (>= " ^ (string_of_int (n + 1)) ^
            " expected)")
    | Some a -> cont a
  with
  | Invalid_argument _ ->
     fail loc "app_to_nth: index must be non-negative"

let app_to_nth_dflt (loc : location) (l : 'a list) (n : int)
  (cont : 'a -> 'b list) : 'b list =
  try app_to_nth loc l n cont with
  | Failure s ->
     print_info loc "%s\n" s;
     []

(* extend current explicit paths with a direction *)
let add_dir (d : dir) (dll : path list) : path list =
  List.map (fun dl -> d :: dl) dll

(* compare literals *)
let is_equal_lit (l : lit) (l' : lit) =
  match l, l' with
  | Lit_unit, Lit_unit -> true
  | Lit_uninitialized, Lit_uninitialized -> true
  | Lit_bool b, Lit_bool b' when b = b' -> true
  | Lit_int n, Lit_int n' when n = n' -> true
  | Lit_double d, Lit_double d' when d = d' -> true
  | Lit_string s, Lit_string s' when s = s' -> true
  | _ -> false

let match_regexp (r : rexp) (t : trm) : bool =
  let aux (r : rexp) (t : trm) : bool =
    let ts = ast_to_string t in
    (* For debug: print on stdout "Considered: %s\n" ts *)
    if r.exact then Str.string_match r.exp ts 0
    else
      try let _ = Str.search_forward r.exp ts 0 in true with
      | Not_found -> false
  in
  if r.only_instr then
    match t.desc with
    | Trm_decl (Def_var _)
      | Trm_apps _
      | Trm_abort _ ->
       t.is_instr && aux r t
    | _ -> false
  else aux r t

(* check if constraint c is satisfied by trm t *)
let rec check_constraint (c : constr) (t : trm) : bool =
  match t.annot with
  | Some Heap_allocated | Some Delete_instructions ->
     (* if t is one of the heap allocation patterns, we simplify it before *)
     check_constraint c (forget_heap_alloc t)
  | Some Access ->
     (* forget the star operator at the root before checking the constraint *)
     begin match t.desc with
     | Trm_apps ({desc = Trm_val (Val_prim (Prim_unop Unop_get)); _}, [t']) ->
        check_constraint c t'
     | _ -> fail t.loc "check_constraint: bad access annotation"
     end
  | Some Multi_decl ->
     (*
       check the constraint on each element of the seq and return true if one
       is true
      *)
     begin match t.desc with
     | Trm_seq tl -> List.mem true (List.map (check_constraint c) tl)
     | _ -> fail t.loc "check_constraint: bad multi_decl annotation"
     end
  | _ ->
     let loc = t.loc in
     begin match c, t.desc with
     (*
       target constraints never hold since they are checked against nodes before
       calling check_constraint in resolve_target
      *)
     | Constr_dir _, _
       | Constr_list _, _
       | Constr_include _, _ ->
        false
     | Constr_regexp r, _ -> match_regexp r t
     | Constr_for (p_init, p_cond, p_step, p_body),
       Trm_for (init, cond, step, body) ->
        check_target p_init init &&
        check_target p_cond cond &&
        check_target p_step step &&
        check_target p_body body
     | Constr_while (p_cond, p_body), Trm_while (cond, body) ->
        check_target p_cond cond &&
        check_target p_body body
     | Constr_if (p_cond, p_then, p_else), Trm_if (cond, then_t, else_t) ->
        check_target p_cond cond &&
        check_target p_then then_t &&
        check_target p_else else_t
     | Constr_decl_var (name, p_body), Trm_decl (Def_var ((x, _), body)) ->
        check_name name x &&
        check_target p_body body
     | Constr_decl_fun (name, cl_args, p_body),
       Trm_decl (Def_fun (x, _, args, body)) ->
        let tl = List.map (fun (x, _) -> trm_var ~loc x) args in
        check_name name x &&
        check_list cl_args tl &&
        check_target p_body body
     | Constr_decl_type name, Trm_decl (Def_typ (x, _)) ->
        check_name name x
     | Constr_decl_enum (name, cec), Trm_decl (Def_enum (n, xto_l)) ->
        check_name name n &&
        check_enum_const cec xto_l
     | Constr_seq cl, Trm_seq tl -> check_list cl tl
     | Constr_var name, Trm_var x -> check_name name x
     | Constr_lit l, Trm_val (Val_lit l') -> is_equal_lit l l'
     | Constr_app (p_fun, cl_args), Trm_apps (f, args) ->
        check_target p_fun f &&
        check_list cl_args args
     | Constr_label (so, p_body), Trm_labelled (l, body) ->
        check_name so l &&
        check_target p_body body
     | Constr_goto so, Trm_goto l -> check_name so l
     | Constr_return p_res, Trm_abort (Ret (Some res)) ->
        check_target p_res res
     | Constr_abort Any, Trm_abort _ -> true
     | Constr_abort Return, Trm_abort (Ret _) -> true
     | Constr_abort Break, Trm_abort Break -> true
     | Constr_abort Continue, Trm_abort Continue -> true
     | Constr_access (p_base, ca), _ ->
        let (base, al) = compute_accesses t in
        check_target p_base base &&
        check_accesses ca al
     | Constr_switch (p_cond, cc), Trm_switch (cond, cases) ->
        check_target p_cond cond &&
        check_cases cc cases
     | _ -> false
     end

and check_name (name : constr_name) (s : string) : bool =
  match name with
  | None -> true
  | Some r -> match_regexp r (trm_var s)

and check_list (cl : constr_list) (tl : trm list) : bool =
  let (p, validate) = cl in
  validate (List.map (check_target p) tl)

and check_accesses (ca : constr_accesses) (al : trm_access list) : bool =
  let rec aux (cal : constr_access list) (al : trm_access list) : bool =
    match cal, al with
    | [], [] -> true
    | Array_access p_index :: cal, Array_access index :: al ->
       check_target p_index index &&
       aux cal al
    | Struct_access so :: cal, Struct_access f :: al ->
       check_name so f &&
       aux cal al
    | Any_access :: cal, _ :: al -> aux cal al
    | _ -> false
  in
  match ca with
  | None -> true
  | Some cal -> aux cal al

and check_cases (cc : constr_cases) (cases : (trm list * trm) list) : bool =
  let rec aux (cl : (case_kind * target) list)
    (cases : (trm list * trm) list) : bool =
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

and check_kind (k : case_kind) (tl : trm list) : bool =
  match k, tl with
  | Case_any, _ -> true
  | Case_default, [] -> true
  | Case_val p_val, _ ->
     (* check that one of the cases corresponds to the given target *)
     List.mem true (List.map (check_target p_val) tl)
  | _ -> false

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

(* check if target tr leads to at least one subterm of t *)
and check_target (tr : target) (t : trm) : bool =
  match resolve_target p t with
  | [] -> false
  | _ -> true

(*
  resolve_target computes the directions to matching subterms
  expected invariants: no duplicate target and no target which is prefix of
  another target that appears after it in the list. Guaranteed by the call to
  sort_unique
 *)
and resolve_target (tr : target)
  (t : trm) : path list =
  let epl =
    match tr with
    | [] -> [[]]
    
    | c :: tr ->
       let dll = resolve_constraint c tr t in
        (explore_in_depth (c :: tr) t) ++ dll
  in
  List.sort_uniq compare_path epl

(* check c against t and in case of success continue with p *)
and resolve_constraint (c : constr) (p : target) (t : trm) : path list =
  let loc = t.loc in
  match c with
  (*
    do not resolve in included files, except if the constraint is Constr_include
   *)
  | Constr_include h when t.annot = Some (Include h) ->
     (*
       remove the include annotation for target resolution to proceed in the
       included file
      *)
     resolve_target p {t with annot = None}
  | _ when is_included t ->
     print_info loc "resolve_constraint: not an include constraint\n";
     []
  (* target constraints first *)
  (* following directions *)
  | Constr_dir d -> follow_dir d p t
  (* list constraint: explore each continuation *)
  | Constr_list (p_elt, next) ->
     (* take into account heap allocation patterns *)
     begin match t.annot with
     | Some Delete_instructions ->
        begin match t.desc with
        (* delete instructions + abort *)
        | Trm_seq ({annot = Some Heap_allocated; _} :: _) ->
           print_info loc
             "resolve_constraint: list constraint applied to a wrong term\n";
           []
        (* instruction + delete instructions *)
        | Trm_seq (t' :: _) -> add_dir (Dir_nth 0) (resolve_constraint c p t')
        | _ -> fail loc "resolve_constraint: bad delete instructions"
        end
     | Some Heap_allocated ->
        (* t is either a heap allocation or a dereferencing *)
        print_info loc
          "resolve_constraint: list constraint applied to a wrong term\n";
        []
     | _ ->
        begin match t.desc with
        | Trm_seq tl ->
           let il = next (List.map (check_target p_elt) tl) in
           explore_list_ind tl (fun n -> Dir_nth n) il (resolve_target p)
        | Trm_apps (_, tl) ->
           let il = next (List.map (check_target p_elt) tl) in
           explore_list_ind tl (fun n -> Dir_arg n) il (resolve_target p)
        | Trm_decl (Def_fun (_, _, args, _)) ->
           let tl = List.map (fun (x, _) -> trm_var ~loc x) args in
           let il = next (List.map (check_target p_elt) tl) in
           explore_list_ind tl (fun n -> Dir_arg n) il (resolve_target p)
        | _ ->
           print_info loc
             "resolve_constraint: list constraint applied to a wrong term\n";
           []
        end
     end
  (*
    if the constraint is a target constraint that does not match the node or
    if it is another kind of constraint, then we check if it holds
   *)
  | c when check_constraint c t -> resolve_target p t
  | _ ->
     print_info loc "resolve_constraint: constraint %s does not hold\n"
       (constraint_to_string c);
     []

(* call resolve_target on subterms of t if possible *)
and explore_in_depth (p : target) (t : trm) : path list =
  let loc = t.loc in
  match t.annot with
  (* no exploration in depth in included files *)
  | Some (Include _) ->
     print_info loc "explore_in_depth: no exploration in included files\n";
     []
  (* we first deal with heap allocation patterns *)
  | Some Delete_instructions ->
     begin match t.desc with
     (* delete instructions + abort *)
     | Trm_seq ({annot = Some Heap_allocated; _} :: tl) ->
        let (n, t') = last tl in
        add_dir (Dir_nth (n + 1)) (explore_in_depth p t')
     (* instruction + delete instructions *)
     | Trm_seq (t' :: _) -> add_dir (Dir_nth 0) (explore_in_depth p t')
     | _ -> fail loc "explore_in_depth: bad delete instructions"
     end
  | Some Heap_allocated ->
     begin match t.desc with
     (* dereferencing *)
     | Trm_apps (_, [t']) -> add_dir (Dir_arg 0) (resolve_target p t')
     (* declaration *)
     | Trm_seq tl ->
        begin match tl with
        (* no initial value (init = uninitialized) *)
        | [{desc = Trm_decl (Def_var ((x, {ty_desc = Typ_ptr _; _}), _)); _}] ->
           add_dir (Dir_nth 0)
             (add_dir Dir_name (resolve_target p (trm_var ~loc x)))
        (* initialisation *)
        | [{desc = Trm_decl (Def_var ((x, _), _)); _} ;
           {desc = Trm_apps (_, [{desc = Trm_var y; _}; init]); _}]
             when y = x ->
           (add_dir (Dir_nth 0)
             (add_dir Dir_name (resolve_target p (trm_var ~ loc x)))) ++
           add_dir (Dir_nth 1) (add_dir (Dir_arg 1) (resolve_target p init))
        | _ -> fail loc "explore_in_depth: bad heap allocation"
        end
     | _ -> fail loc "explore_in_depth: bad heap_alloc instruction"
     end
  | Some Access ->
     begin match t.desc with
       (*
         the wildcard is a star operator the user doesn't know about
         t' is an access under which want to explore
        *)
     | Trm_apps (_, [t']) -> add_dir (Dir_arg 0) (explore_in_depth p t')
     | _ -> fail loc "explore_in_depth: bad access annotation"
     end
  | Some Multi_decl ->
     (* explore each declaration in the seq *)
     begin match t.desc with
     | Trm_seq tl ->
        explore_list tl (fun i -> Dir_nth i) (explore_in_depth p)
     | _ -> fail loc "explore_in_depth: bad multi_decl annotation"
     end
  | _ ->
     begin match t.desc with
     | Trm_decl (Def_var ((x, _), body))
       | Trm_decl (Def_fun (x, _, _, body)) ->
        add_dir Dir_name (resolve_target p (trm_var ~loc x)) ++
        add_dir Dir_body (resolve_target p body)
     | Trm_decl (Def_enum (x, xto_l)) ->
        let (il, tl) =
          foldi
            (fun n (il, tl) (_, t_o) ->
              match t_o with
              | None -> (il, tl)
              | Some t -> (il ++ [n], tl ++ [t])
            )
           ([], [])
           xto_l
        in
        add_dir Dir_name (resolve_target p (trm_var ~loc x)) ++
        (explore_list (List.map (fun (y, _) -> trm_var ~loc y) xto_l)
           (fun n -> Dir_enum_const (n, Enum_const_name))
           (resolve_target p)) ++
        (explore_list tl
           (fun n -> Dir_enum_const (List.nth il n, Enum_const_val))
           (resolve_target p))
     | Trm_abort (Ret (Some body)) ->
        add_dir Dir_body (resolve_target p body)
     | Trm_for (init, cond, step, body) ->
        (* init *)
        (add_dir Dir_for_init (resolve_target p init)) ++
        (* cond *)
        (add_dir Dir_cond (resolve_target p cond)) ++
        (* step *)
        (add_dir Dir_for_step (resolve_target p step)) ++
        (* body *)
        (add_dir Dir_body (resolve_target p body))
     | Trm_while (cond, body) ->
        (* cond *)
        (add_dir Dir_cond (resolve_target p cond)) ++
        (* body *)
        (add_dir Dir_body (resolve_target p body))
     | Trm_if (cond, then_t, else_t) ->
        (* cond *)
        (add_dir Dir_cond (resolve_target p cond)) ++
        (* then *)
        (add_dir Dir_then (resolve_target p then_t)) ++
        (* else *)
        (add_dir Dir_else (resolve_target p else_t))
     | Trm_apps (f, args) ->
        (* fun *)
        (add_dir Dir_app_fun (resolve_target p f)) ++
        (* args *)
        (explore_list args (fun n -> Dir_arg n) (resolve_target p))
     | Trm_seq tl
       | Trm_array tl
       | Trm_struct tl ->
        explore_list tl (fun n -> Dir_nth n) (resolve_target p)
     | Trm_val (Val_array vl)
       | Trm_val (Val_struct vl) ->
        explore_list (List.map (trm_val ~loc) vl) (fun n -> Dir_nth n)
          (resolve_target p)
     | Trm_labelled (l, body) ->
        add_dir Dir_name (resolve_target p (trm_var ~loc l)) ++
        add_dir Dir_body (resolve_target p body)
     | Trm_switch (cond, cases) ->
        (add_dir Dir_cond (resolve_target p cond)) ++
        (foldi (fun i epl case -> epl ++ explore_case i case p) [] cases)
     | _ ->
        print_info loc "explore_in_depth: cannot find a subterm to explore\n";
        []
     end

(*
  call resolve_target on given case name and body
  i is the index of the case in its switch
 *)
and explore_case (i : int) (case : trm list * trm) (p : target) : path list =
  let (tl, body) = case in
  match tl with
  (* default case *)
  | [] ->
     add_dir (Dir_case (i, Case_body)) (resolve_target p body)
  | _ ->
     (foldi
        (fun j epl t ->
          epl ++
          (add_dir (Dir_case (i, Case_name j)) (resolve_target p t))
        )
        []
        tl
     ) ++
     add_dir (Dir_case (i, Case_body)) (resolve_target p body)

(* follow the direction d in t and call resolve_target on p *)
and follow_dir (d : dir) (p : target) (t : trm) : path list =
  let loc = t.loc in
  match d, t.desc with
  | Dir_nth n, Trm_seq tl
    | Dir_nth n, Trm_array tl
    | Dir_nth n, Trm_struct tl ->
     app_to_nth_dflt loc tl n
       (fun nth_t -> add_dir (Dir_nth n) (resolve_target p nth_t))
  | Dir_nth n, Trm_val (Val_array vl)
    | Dir_nth n, Trm_val (Val_struct vl) ->
     app_to_nth_dflt loc vl n (fun nth_v ->
         add_dir (Dir_nth n) (resolve_target p (trm_val ~loc nth_v)))
  | Dir_cond, Trm_if (cond, _, _)
    | Dir_cond, Trm_while (cond, _)
    | Dir_cond, Trm_for (_, cond, _, _)
    | Dir_cond, Trm_switch (cond, _) ->
     add_dir Dir_cond (resolve_target p cond)
  | Dir_then, Trm_if (_, then_t, _) ->
     add_dir Dir_then (resolve_target p then_t)
  | Dir_else, Trm_if (_, _, else_t) ->
     add_dir Dir_else (resolve_target p else_t)
  | Dir_body, Trm_decl (Def_var (_, body))
    | Dir_body, Trm_decl (Def_fun (_, _, _, body))
    | Dir_body, Trm_for (_, _, _, body)
    | Dir_body, Trm_while (_, body)
    | Dir_body, Trm_abort (Ret (Some body))
    | Dir_body, Trm_labelled (_, body) ->
     add_dir Dir_body (resolve_target p body)
  | Dir_for_init, Trm_for (init, _, _, _) ->
     add_dir Dir_for_init (resolve_target p init)
  | Dir_for_step, Trm_for (_, _, step, _) ->
     add_dir Dir_for_step (resolve_target p step)
  | Dir_app_fun, Trm_apps (f, _) -> add_dir Dir_app_fun (resolve_target p f)
  | Dir_arg n, Trm_apps (_, tl) ->
     app_to_nth_dflt loc tl n (fun nth_t ->
         add_dir (Dir_arg n) (resolve_target p nth_t))
  | Dir_arg n, Trm_decl (Def_fun (_, _, arg, _)) ->
     let tl = List.map (fun (x, _) -> trm_var ~loc x) arg in
     app_to_nth_dflt loc tl n (fun nth_t ->
         add_dir (Dir_arg n) (resolve_target p nth_t))
  | Dir_name, Trm_decl (Def_var ((x, _), _))
    | Dir_name, Trm_decl (Def_fun (x, _, _, _))
    | Dir_name, Trm_decl (Def_enum (x, _))
    | Dir_name, Trm_labelled (x, _)
    | Dir_name, Trm_goto x ->
     add_dir Dir_name (resolve_target p (trm_var ~loc x))
  | Dir_case (n, cd), Trm_switch (_, cases) ->
     app_to_nth_dflt loc cases n
       (fun (tl, body) ->
         match cd with
         | Case_body ->
            add_dir (Dir_case (n, cd)) (resolve_target p body)
         | Case_name i ->
            app_to_nth_dflt loc tl i (fun ith_t ->
                add_dir (Dir_case (n, cd)) (resolve_target p ith_t))
       )
  | Dir_enum_const (n, ecd), Trm_decl (Def_enum (_, xto_l)) ->
     app_to_nth_dflt loc xto_l n
       (fun (x, t_o) ->
         match ecd with
         | Enum_const_name ->
            add_dir (Dir_enum_const (n, ecd)) (resolve_target p (trm_var ~loc x))
         | Enum_const_val ->
            begin match t_o with
            | None ->
               print_info loc "follow_dir: no value for constant of index %d\n"
                 n;
               []
            | Some t ->
               add_dir (Dir_enum_const (n, ecd)) (resolve_target p t)
            end
       )
  | _, _ ->
     print_info loc "follow_dir: direction %s does not match"
       (dir_to_string d);
     []

(*
  call cont on each element of the list and gathers the results
  used for seq, array, struct, and fun arguments
  d: function that gives the direction to add depending on the index
 *)
and explore_list (tl : trm list) (d : int -> dir)
  (cont : trm -> path list) : path list =
  foldi (fun i epl t -> epl ++ add_dir (d i) (cont t)) [] tl

(*
  call cont on each element of the list whose index is in the domain and
  gather the results
  d: function that gives the direction to add depending on the index
 *)
and explore_list_ind (tl : trm list) (d : int -> dir) (dom : int list)
  (cont : trm -> path list) : path list =
  foldi
    (fun i epl t ->
      if List.mem i dom then epl ++ add_dir (d i) (cont t) else epl)
    []
    tl

(******************************************************************************)
(*                         Explicit target manipulation                         *)
(******************************************************************************)

(*
  follow the explicit path and return the corresponding subterm and its context
 *)
let resolve_explicit_path (dl : path) (t : trm) : trm * (trm list) =
  let rec aux (dl : path) (t : trm) (ctx : trm list) : trm * (trm list) =
    match dl with
    | [] -> (t, List.rev ctx)
    | d :: dl ->
       let loc = t.loc in
       begin match d, t.desc with
       | Dir_nth n, Trm_seq tl ->
          let decl_before (n : int) (tl : trm list) =
            foldi
              (fun i acc (t : trm) ->
                if i >= n then acc
                else
                  match t.desc with
                  | Trm_decl _ -> t :: acc
                  | Trm_seq _ when t.annot = Some Heap_allocated -> t :: acc
                  | _ -> acc
              )
              []
              tl
          in
          app_to_nth loc tl n
            (fun nth_t -> aux dl nth_t ((decl_before n tl) ++ ctx))
       | Dir_nth n, Trm_array tl
         | Dir_nth n, Trm_struct tl ->
          app_to_nth loc tl n (fun nth_t -> aux dl nth_t ctx)
       | Dir_nth n, Trm_val (Val_array vl)
         | Dir_nth n, Trm_val (Val_struct vl) ->
          app_to_nth loc vl n
            (fun v -> aux dl (trm_val ~loc v) ctx)
       | Dir_cond, Trm_if (cond, _, _)
         | Dir_cond, Trm_while (cond, _)
         | Dir_cond, Trm_switch (cond, _) ->
          aux dl cond ctx
       | Dir_cond, Trm_for (init, cond, _, _) ->
          begin match init.desc with
          (* loop indices are heap allocated *)
          | Trm_seq _ when init.annot = Some Heap_allocated ->
             aux dl cond (init :: ctx)
          | _ -> aux dl cond ctx
          end
       | Dir_then, Trm_if (_, then_t, _) ->
          aux dl then_t ctx
       | Dir_else, Trm_if (_, _, else_t) ->
          aux dl else_t ctx
       | Dir_body, Trm_decl (Def_fun (_, _, args, body)) ->
          (* do as if fun args were heap allocated *)
          let args_decl =
            List.rev_map
              (fun (x, tx) ->
                trm_seq ~annot:(Some Heap_allocated)
                  [trm_decl (Def_var ((x, typ_ptr tx),
                                      trm_lit Lit_uninitialized))]
              )
              args
          in
          aux dl body (args_decl ++ ctx)
       | Dir_body, Trm_for (init, _, _, body) ->
          begin match init.desc with
          | Trm_seq _ when init.annot = Some Heap_allocated ->
             aux dl body (init :: ctx)
          | _ -> aux dl body ctx
          end
       | Dir_body, Trm_decl (Def_var (_, body))
         | Dir_body, Trm_while (_, body)
         | Dir_body, Trm_abort (Ret (Some body))
         | Dir_body, Trm_labelled (_, body) ->
          aux dl body ctx
       | Dir_for_init, Trm_for (init, _, _, _) ->
          aux dl init ctx
       | Dir_for_step, Trm_for (init, _, step, _) ->
          begin match init.desc with
          | Trm_seq _ when init.annot = Some Heap_allocated ->
             aux dl step (init :: ctx)
          | _ -> aux dl step ctx
          end
       | Dir_app_fun, Trm_apps (f, _) -> aux dl f ctx
       | Dir_arg n, Trm_apps (_, tl) ->
          app_to_nth loc tl n (fun nth_t -> aux dl nth_t ctx)
       | Dir_arg n, Trm_decl (Def_fun (_, _, arg, _)) ->
          app_to_nth loc arg n
            (fun (x, _) -> aux dl (trm_var ~loc x) ctx)
       | Dir_name, Trm_decl (Def_var ((x, _), _))
         | Dir_name, Trm_decl (Def_fun (x, _, _, _))
         | Dir_name, Trm_decl (Def_enum (x, _))
         | Dir_name, Trm_labelled (x, _)
         | Dir_name, Trm_goto x ->
          aux dl (trm_var ~loc x) ctx
       | Dir_case (n, cd), Trm_switch (_, cases) ->
          app_to_nth loc cases n
            (fun (tl, body) ->
              match cd with
              | Case_body -> aux dl body ctx
              | Case_name i ->
                 app_to_nth loc tl i (fun ith_t -> aux dl ith_t ctx)
            )
       | Dir_enum_const (n, ecd), Trm_decl (Def_enum (_, xto_l)) ->
          app_to_nth loc xto_l n
             (fun (x, t_o) ->
               match ecd with
               | Enum_const_name -> aux dl (trm_var ~loc x) ctx
               | Enum_const_val ->
                  begin match t_o with
                  | None ->
                     fail loc
                       "resolve_explicit_path: no value for enum constant"
                  | Some t ->
                     aux dl t ctx
                  end
             )
       | _, _ ->
          let s = dir_to_string d in
          fail loc ("resolve_explicit_path: direction " ^ s ^ " does not match")
       end
  in
  aux dl t []

(*
  find the explicit path to the toplevel declaration of x if it exists
  assumption: x denotes a function or a type
  todo: generalise to other terms
 *)
let rec path_to_decl (x : var) (t : trm) : path option =
  match t.desc with
  | Trm_decl def ->
     begin match def with
     | Def_fun (f, _, _, _) when f = x -> Some []
     | Def_typ (y, _) when y = x -> Some []
     | _ -> None
     end
  | Trm_seq tl ->
     foldi
       (fun i dlo t' ->
         match dlo with
         | Some _ -> dlo
         | _ ->
            begin match path_to_decl x t' with
            | Some dl -> Some (Dir_nth i :: dl)
            | _ -> None
            end
       )
       None
       tl
  (* val, var, array, struct, if, apps, while, for, switch, abort, label *)
  | _ -> None
