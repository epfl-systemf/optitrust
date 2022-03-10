open Ast

(* [intro_mcalloc tg] expects the target [tg] pointing to a call to funciton alloc
      then it will replace this call with a call to MCALLOC
*)
let intro_mcalloc : Target.Transfo.t =
  Target.apply_on_targets (Matrix_core.intro_mcalloc)

(* [intro_mmalloc tg] expects the target [tg] pointing to a call to funciton malloc
      then it will replace this call with a call to MMALLOC
*)
let intro_mmalloc : Target.Transfo.t =
  Target.apply_on_targets (Matrix_core.intro_mmalloc)

(* [intro_mindex dim tg] expects the target [tg] pointing to an array access
    then it will replace that access to let say index i with an access at
    MINDEX (dim,i)
*)
let intro_mindex (dim : trm) : Target.Transfo.t =
  Target.apply_on_targets (Matrix_core.intro_mindex dim)

(* [reorder_dims order tg]: expects the target [tg] pointing to a call to ALLOC functions, or MINDEX
      then it will reorder their args based on [order], where [order] is a list of indices which the
      current args should follow
*)
let reorder_dims ?(rotate_n : int = 0) ?(order : int list = [])  (): Target.Transfo.t =
  Target.apply_on_targets (Matrix_core.reorder_dims rotate_n order)

(* [insert_alloc_dim new_dim]: expects the target [tg] pointing to call to ALLOC functions, then it will
      add a new arg at the begining of the list of args in the targeted call
 *)
let insert_alloc_dim (new_dim : trm) : Target.Transfo.t =
  Target.apply_on_targets (Matrix_core.insert_alloc_dim new_dim)

(* [insert_access_dim new_dim new_index tg] expects the target [tg] pointing to an array access, then it will
    add two new args in the call to MINDEX function inside that array access *)

let insert_access_dim_index (new_dim : trm) (new_index : trm) : Target.Transfo.t =
  Target.apply_on_targets (Matrix_core.insert_access_dim_index new_dim new_index)

(* [biject fun_name tg] expectes the target [tg] poiting to a function call then it
      replaces the name of that function call with [fun_name]
*)
let biject (fun_name : string) : Target.Transfo.t =
  Expr.replace_fun fun_name

(* [local_name ~mark var into tg] expects the target pointing to an instruction that contains
      an occurrence of [var] then it will define a matrix [into] whose dimensions will be the same
      as the one of [var]. Then we copy the contents of the matrix [var] into [into] and finaly we
      free up the memory.
 *)

let local_name ?(my_mark : mark option) ?(indices : (var list) = []) ?(alloc_target : Target.target option = None) (v : var) ~into:(into : var) (tg : Target.target) : unit =
  let remove = (my_mark = None) in 
  Internal.nobrace_remove_after ~remove (fun _ -> 
    Target.(apply_on_targets 
     (fun t p -> 
        let seq_p, i = Internal.isolate_last_dir_in_seq p in
        let seq = Target.target_of_path seq_p in 
        let var_target = cOr [[cVarDef v];[cWriteVar v]] in 
        let vardef_trm = Target.get_trm_at (seq @ [var_target]) in 
        let var_type = match vardef_trm.desc with 
          | Trm_let (_, (_, ty), _) -> get_inner_ptr_type ty 
          | Trm_apps (_, [lhs; _rhs]) when is_set_operation vardef_trm ->
            begin match lhs.typ with 
            | Some ty -> ty
            | None -> fail vardef_trm.loc "local_name: couldn't find the type of the targetd variable'"
            end
          | _ ->fail vardef_trm.loc "local_name: couldn't find the type of the targetd variable'"
        in 

        let alloc_trm = if alloc_target <> None then Target.get_trm_at (Tools.unsome alloc_target)
          else Target.(get_trm_at (seq @ [var_target; Target.cFun ~regexp:true "M.ALLOC."])) in 
        
        let alloc_trms = match Matrix_core.alloc_inv alloc_trm with
        | Some (dims, sz, zero_init) -> (dims, sz, zero_init)
        | _ -> fail None "local_name: could not get the dimensions and the size of the matrix" in
         
        if not remove then Internal.nobrace_enter();
        Matrix_core.local_name my_mark v into alloc_trms var_type indices t p
     ) tg)
  ) 



let local_name1 ?(my_mark : mark option) ?(indices : (var list) = []) ?(is_detached : bool = false) (var : var) ~into:(into : var) (tg : Target.target) : unit =
  let vardef_trm = Target.get_trm_at [Target.cVarDef var] in
  let var_type = match trm_var_def_inv vardef_trm with
  | Some (_, _, ty, _) -> ty
  | _ -> fail vardef_trm.loc "local_name: make sure the name of the current var is entered correctly" in
  let alloc_tg = if not is_detached then  [Target.cVarDef var; Target.cFun ~regexp:true "M.ALLOC."] else [Target.cWriteVar var; Target.cFun ~regexp:true "M.ALLOC."] in 
  let alloc_trm = Target.get_trm_at alloc_tg in 
  (* let alloc_trm = Target.get_trm_at [Target.cVarDef var; Target.cFun ~regexp:true "M.ALLOC."] in *)
  let alloc_trms = match Matrix_core.alloc_inv alloc_trm with
  | Some (dims, sz, zero_init) -> (dims, sz, zero_init)
  | _ -> fail None "local_name: could not get the dimensions and the size of the matrix" in
  begin match my_mark with
  | Some _ -> Internal.nobrace_enter (); Target.apply_on_targets (Matrix_core.local_name my_mark var into alloc_trms var_type indices) tg
  | _ ->
  Internal.nobrace_remove_after (fun _ ->
    Target.apply_on_targets (Matrix_core.local_name my_mark var into alloc_trms var_type indices ) tg
  ) end


(* [delocalize ~init_zero ~acc_in_place ~acc ~dim ~index ~ops] a generalized version of variable_delocalize*)
let delocalize ?(init_zero : bool = false) ?(acc_in_place : bool = false) ?(acc : string option) ?(any_mark : mark = "")~dim:(dim : trm)  ~index:(index : string) ~ops:(dl_o : delocalize_ops) : Target.Transfo.t =
    Target.apply_on_targets (Matrix_core.delocalize dim init_zero acc_in_place acc any_mark index dl_o)