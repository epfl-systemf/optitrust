open Ast
open Target
open Clang_to_ast
open Ast_to_c
open Tools
(* INTERNAL FUNCTIONS *)
(* *********************************************** *)
(*replace occurrences of t_before with t_after in t
  paths point at subterms in which all occurences will be replaced
  the empty path means all occurences will be replaced (default behaviour)
  assumption: t_before and t_after are equivalent (in terms of value and of side
  effects)
 *)
let change_trm ?(change_at : target list = [[]]) (t_before : trm)
  (t_after : trm) (t : trm) : trm =
  (* change all occurences of t_before in t' *)
  let rec apply_change (t' : trm) =
    (* necessary because of annotations that may be different *)
    if ast_to_string t' = ast_to_string t_before then t_after
    else
      match t'.desc with
      (*
        particular case for heap allocation: do not change the lhs of the
        initialisation
       *)
      | Trm_seq [t_decl; {desc = Trm_apps (_, [lhs; init]); loc; _}]
           (* when t'.annot = Some Heap_allocated *) ->
         trm_seq ~annot:t'.annot ~loc:t'.loc ~add:t'.add
           ~attributes:t'.attributes
           [
             t_decl;
             trm_set (* ~annot:(Some Initialisation_instruction) *) ~loc lhs
               (apply_change init)
           ]
      | _ -> trm_map apply_change t'
  in
  List.fold_left
    (fun t' tr ->
      let b = !Flags.verbose in
      Flags.verbose := false;
      let epl = resolve_target tr t' in
      Flags.verbose := b;
      match epl with
      | [] ->
         print_info t'.loc "change_trm: no matching subterm for target %s\n"
           (target_to_string tr);
         t'
      | _ -> List.fold_left (apply_on_path apply_change) t' epl
    )
    t
    change_at



(* same as change_trm but for types *)
let change_typ ?(change_at : target list = [[]]) (ty_before : typ)
  (ty_after : typ) (t : trm) : trm =
  (* change all occurences of ty_before in ty *)
  let rec change_typ (ty : typ) : typ =
    (* necessary because of annotations in trms that may be different *)
    if typ_to_string ty = typ_to_string ty_before then ty_after
    else Ast.typ_map change_typ ty
  in
  (* change all occurrences of ty_before in type annotations in t *)
  let rec replace_type_annot (t : trm) : trm =
    let t =
      {t with typ = match t.typ with
                    | None -> None
                    | Some ty' -> Some (change_typ ty')
      }
    in
    trm_map replace_type_annot t
  in
  (* change all occurences of ty_before in t *)
  let apply_change (t : trm) : trm =
    let rec aux (t : trm) : trm =
      (* only match nodes where typs occur *)
      match t.desc with
      | Trm_val (Val_prim (Prim_new ty)) ->
         trm_prim ~annot:t.annot ~loc:t.loc ~add:t.add
           (Prim_new (change_typ ty))
      | Trm_val (Val_prim (Prim_unop (Unop_cast ty))) ->
         trm_unop ~annot:t.annot ~loc:t.loc ~add:t.add
           (Unop_cast (change_typ ty))
      | Trm_let (vk,(y,ty),init) ->
        trm_let ~annot:t.annot ~loc:t.loc ~is_statement:t.is_statement ~add:t.add vk (y,change_typ ty) (aux init)
      | Trm_let_fun (f, ty, args, body) ->
         trm_let_fun ~annot:t.annot ~loc:t.loc ~is_statement:t.is_statement ~add:t.add
           ~attributes:t.attributes f (change_typ ty)
                     (List.map (fun (y, ty) -> (y, change_typ ty)) args)
                     (aux body)
      | Trm_typedef (Typedef_abbrev (y, ty)) ->
         trm_typedef  ~annot:t.annot ~loc:t.loc ~is_statement:t.is_statement ~add:t.add ~attributes:t.attributes 
          (Typedef_abbrev (y, change_typ ty))
      | _ -> trm_map aux t
    in
    replace_type_annot (aux t)
  in
  List.fold_left
    (fun t' tr ->
      let b = !Flags.verbose in
      Flags.verbose := false;
      let epl = resolve_target tr t' in
      Flags.verbose := b;
      match epl with
      | [] ->
         print_info t'.loc "change_typ: no matching subterm for target %s\n"
           (target_to_string tr);
         t'
      | _ -> List.fold_left (apply_on_path apply_change) t' epl
    )
    t

    change_at

(*
  insert inert after the subterm pointed at by dl in t
  assumption: dl points at a seq element, thus ends with Dir_nth n
  if the inserted element must be first in the seq, use n < 0
 *)
let insert_trm_after (dl : path) (insert : trm) (t : trm) : trm =
  let dl' = List.rev dl in
  match List.hd dl' with
  | Dir_nth n ->
     apply_on_path
       (fun t' ->
         match t'.desc with
         | Trm_seq tl ->
            trm_seq ~annot:t'.annot ~loc:t'.loc ~add:t'.add
              ~attributes:t'.attributes (Tools.list_insert n insert tl)
         | _ -> fail t'.loc "insert_trm_after: path points at wrong term"
       )
       t
       (List.rev (List.tl dl'))
  | _ -> fail t.loc "insert_trm_after: bad path"

(* make sure each occurence of y in t is marked with type variable x *)
let rec replace_type_with (x : typvar) (y : var) (t : trm) : trm =
  match t.desc with
  | Trm_var y' when y' = y ->
     trm_var ~annot:t.annot ~loc:t.loc ~add:t.add ~typ:(Some (typ_var x (get_typedef x)))
       ~attributes:t.attributes y
  | _ -> trm_map (replace_type_with x y) t



(*
  replace with x the types of the variables given by their index
  assumption: t is a fun body whose arguments are given by tvl
 *)
let replace_arg_types_with (x : typvar) (il : int list) (tvl : typed_var list)
  (t : trm) : trm =
  List.fold_left
    (fun t' i ->
      let (y, _) = List.nth tvl i in
      replace_type_with x y t'
    )
    t
    il


let rec functions_with_arg_type ?(outer_trm : trm option = None) (x : typvar)
  (t : trm) : ilset funmap =
  let rec aux (t : trm) : ilset Tools.funmap =
    match t.desc with
    | Trm_let (_,(_,_),body) -> aux body
    | Trm_let_fun (_, _, _, body) -> aux body
    | Trm_if (cond, then_, else_) -> aux cond +@ aux then_ +@ aux else_
    | Trm_seq tl ->
       List.fold_left (fun ilsm t' -> ilsm +@ aux t') Fun_map.empty tl
    | Trm_apps (f, tl) ->
       (* functions may be applied to arguments of type x in terms of tl *)
       let ilsm =
         List.fold_left (fun ilsm t' -> ilsm +@ aux t') Fun_map.empty tl
       in
       begin match f.desc with
       (*
         if f is a variable, we have to add f to ilsm if an argument has type x
         ignore the free function
        *)
       | Trm_var f when f <> "free" ->
          let il =
            foldi
              (fun i il (t' : trm) ->
                match t'.typ with
                (* note: also works for heap allocated variables *)
                | Some {ty_desc = Typ_var (x', _); _} when x' = x -> i :: il
                | _ -> il
              )
              []
              tl
          in
          begin match il with
          | [] -> ilsm
          | _ ->
             let ils = IntListSet.singleton (List.rev il) in
             Fun_map.update f
               (function
                | None -> Some ils
                | Some ils' -> Some (IntListSet.union ils ils')
               )
               ilsm
          end
       (* in other cases, do a recursive call *)
       | _ -> ilsm +@ aux f
       end
    | Trm_while (cond, body) -> aux cond +@ aux body
    | Trm_for (init, cond, step, body) ->
       aux init +@ aux cond +@ aux step +@ aux body
    | Trm_switch (cond, cases) ->
       aux cond +@
         List.fold_left (fun ilsm t' -> ilsm +@ aux t') Fun_map.empty
           (* no function applications in case values *)
           (List.map (fun (_, t') -> t') cases)
    | Trm_abort (Ret (Some t'))
      | Trm_labelled (_, t') ->
       aux t'
    (* val, var, array, struct, type decl, aborts with no argument *)
    | _ -> Fun_map.empty
  in
  let ilsm = aux t in
  (*
    for each function, do a recursive call on its declaration where the type of
    arguments is replaced with x
   *)
  Fun_map.fold
    (fun f ils res ->
      IntListSet.fold
        (fun il res ->
          (*
            first compute the body of f where the arguments at positions in il
            have type x
           *)
          let global_trm =
            match outer_trm with
            | None -> t
            | Some t' -> t'
          in
          match target_to_decl f global_trm with
          (* if the declaration cannot be found, ignore this function *)
          | None ->
             print_info global_trm.loc
               ("functions_with_arg_type: cannot find declaration of " ^^
                  "function %s, ignoring it.\n") f;
             Fun_map.remove f res
          | Some dl ->
             let (def, _) = resolve_path dl global_trm in
             begin match def.desc with
             | Trm_let_fun (_, _, args, body) ->
                let b = replace_arg_types_with x il args body in
                (* then do a recursive call on the new body *)
                res +@ functions_with_arg_type ~outer_trm:(Some global_trm) x b
             | _ ->
                fail t.loc
                  ("functions_with_arg_type: wrong target to declaration of " ^ f)
             end
        )
        ils
        res
    )
    ilsm
    ilsm

let clean_up_no_brace_seq (t : trm) : trm =
  let rec clean_up_in_list (tl : trm list) : trm list =
    match tl with
    | [] -> []
    | t :: tl ->
       begin match t.desc with
       | Trm_seq tl' when t.annot = Some No_braces ->
          tl' ++ (clean_up_in_list tl)
       | _ -> t :: (clean_up_in_list tl)
       end
  in
  let rec aux (t : trm) : trm =
    match t.desc with
    (*
      condition on annotation: toplevel insdeclarations might contain a heap
      allocated variable and hence we can find a no_brace seq inside a
      delete_instructions seq, which we do not want to inline
     *)
    | Trm_seq tl (* when t.annot <> Some Delete_instructions *) ->
       trm_seq ~annot:t.annot ~loc:t.loc ~add:t.add ~attributes:t.attributes
         (clean_up_in_list (List.map aux tl))
    | _ -> trm_map aux t
  in
  aux t

(*
  add copies of the provided functions given by their name
  each function f is mapped to the set of lists of indices corresponding to uses
  of f where the arguments at those indices have x for type
  the name of the copies are indexed with "new_name_i", where "new_name" is the
  result of name, and similarly for labels in these copies
 *)
let rec insert_fun_copies (name : var -> var) (ilsm : ilset funmap) (x : typvar)
  (t : trm) : trm =
  (* also change the labels in the body of fun copies for uniqueness *)
  let rec label_aux (i : int) (t : trm) : trm =
    match t.desc with
    | Trm_labelled (l, body) ->
       trm_labelled ~annot:t.annot ~loc:t.loc ~add:t.add
         ~attributes:t.attributes (name l ^ "_" ^ string_of_int i)
         (label_aux i body)
    | _ -> trm_map (label_aux i) t
  in
  clean_up_no_brace_seq
    (Fun_map.fold
       (fun f ils t' ->
         match target_to_decl f t' with
         | None ->
            fail t'.loc
              ("insert_fun_copies: cannot find declaration of function " ^ f)
         | Some dl ->
            let (fdecl, _) = resolve_path dl t' in
            begin match fdecl.desc with
            | Trm_let_fun (f', r, tvl, b) when f = f' ->
               (* for each element of ils, create a copy *)
               let tl =
                 intl_set_foldi
                   (fun i il tl ->
                     (*
                       for each argument whose index is in il, use x as
                       (possibly new) type in the declaration
                      *)
                     let tvl' =
                       List.fold_left
                         (change_nth (fun (y, _) -> (y, typ_var x (get_typedef x)))) tvl il
                     in
                     (* add index to labels in the body of the function *)
                     let b' =
                       label_aux i (replace_arg_types_with x il tvl' b)
                     in
                     (* create the copy of f corresponding to il *)
                     trm_let_fun (name f ^ "_" ^ string_of_int i) r tvl' b' :: tl
                   )
                   ils
                   []
               in
               (* insert the copies of f *)
               insert_trm_after dl
                 (trm_seq ~annot:(Some No_braces) (List.rev tl)) t'
            | _ -> fail t'.loc "insert_fun_copies: bad target to fun decl"
            end
       )
       ilsm
       t
    )

(*
  replace the applications of the given functions to arguments of type x with
  the application of their copies whose name is given by name
 *)
and replace_fun_names (name : var -> var) (ilsm : ilset funmap) (x : typvar)
  (t : trm) : trm =
  let annot = t.annot in
  let loc = t.loc in
  let is_statement = t.is_statement in
  let add = t.add in
  let typ = t.typ in
  let attributes = t.attributes in
  match t.desc with
  | Trm_apps (fun_, args) ->
     begin match fun_.desc with
     | Trm_var f ->
        (* first check if f is one of the functions that required copies *)
        begin match Fun_map.find_opt f ilsm with
        | None -> (* if not, just do a recursive call on args *)
           trm_map (replace_fun_names name ilsm x) t
        | Some ils ->
           (*
             if f required copies, compute the indices of arguments of type x
            *)
           let il =
             List.rev
               (foldi
                  (fun i il (ti : trm) ->
                    match ti.typ with
                    (* note: also works for heap allocated variables *)
                    | Some {ty_desc = Typ_var (x', _); _} when x' = x -> i :: il
                    | _ -> il
                  )
                  []
                  args
               )
           in
           begin match il with
           (* if il = [] then no argument is of type x so do a recursive call *)
           | [] -> trm_map (replace_fun_names name ilsm x) t
           | _ -> (* otherwise, find the appropriate name *)
              let io =
                intl_set_foldi
                  (fun i il' io ->
                    match io with
                    | Some _ -> io
                    | None ->
                       if IntList.compare il il' = 0 then Some i else None
                  )
                  ils
                  None
              in
              let f' =
                match io with
                | None -> fail loc "replace_fun_names: unmatched call"
                | Some i -> name f ^ "_" ^ string_of_int i
              in
              (* also do a recursive call on args *)
              let args' = List.map (replace_fun_names name ilsm x) args in
              trm_apps ~annot ~loc ~is_statement ~add ~typ ~attributes
                (trm_var ~annot:fun_.annot ~loc:fun_.loc ~add:fun_.add
                   ~attributes:fun_.attributes f') args'
           end
        end
     | _ -> trm_map (replace_fun_names name ilsm x) t
     end
  | _ -> trm_map (replace_fun_names name ilsm x) t

(* [isolate_last_dir_in_seq dl]:  
    params: 
      dl: explicit path to the targeted trm
    return:
      a pair of the explicit path to the outer sequence and the index of the term inside that sequence
*)
let isolate_last_dir_in_seq (dl : path) : path * int = 
  match List.rev dl with
  | Dir_nth i :: dl' -> (List.rev dl',i)
  | _ -> fail None "isolate_last_dir_in_seq: cannot isolate the definition in a sequence"



(* ********************************************** *)






(* [var_init_detach_aux t]: This is an auxiliary function for var_init_detach
    params:
      t: an ast subterm
    return:
      a sequence which contains the declaration of the variable and a set operations for that variable
*)
let var_init_detach_aux (t : trm) : trm =
  match t.desc with 
  | Trm_let(vk,(x, tx), init) ->
    begin match vk with 
    | Var_immutable -> fail t.loc "var_init_detach_aux: const declarations cannot be detached"
    | _ ->
      let init = 
        begin match init.desc with 
        | Trm_apps(_,[init]) -> init
        | _ -> fail t.loc "var_init_detach_aux: expected a heap allocated variable declaration"
        end in
      trm_seq [
        trm_let vk (x, tx) (trm_prim (Prim_new tx));
        trm_set (trm_var x) init
      ]
    end
  | _ -> fail t.loc "var_init_detach_aux: variable could not be matched, make sure your path is correct"


let var_init_detach : Target.Transfo.local =
  Target.apply_on_path(var_init_detach_aux ) 

(* [var_init_attach_aux t]: This is an auxiliary function for var_init_attach
    params:
      t: an ast subterm
    return
      the updated
*)
let var_init_attach_aux (t : trm) : trm =
  match t.desc with 
  | Trm_seq tl -> 
    (* Assumption: The sequence is of the form
      {
        int x;
        x = 5;
      }
    *)
    let var_decl = List.nth tl 0 in
    let var_set = List.nth tl 1 in
    let var_kind, var_name, var_type = begin match var_decl.desc with 
    | Trm_let (vk,(x,tx),_) -> vk, x, tx
    | _ -> fail t.loc "var_init_attach_aux: sequence does not satisfy the assumption described above"
    end
    in
    let var_init = begin match var_set.desc with 
    | Trm_apps(_, [_;init]) -> init
    | _ -> fail t.loc "var_init_attach_aux: sequence does not satisfy the assumtion that the second term of the sequence is the set operations"
    end
    in
    trm_let ~loc:t.loc var_kind (var_name, var_type) (trm_apps (trm_prim ~loc:t.loc (Prim_new var_type)) [var_init])
  | _ -> fail t.loc "var_init_attach_aux: sequence was not matched, make sure the path is correct"

(* [var_init_attach t]: Change a sequence of the form {int x; x = 5;} to int x = 5 
    params:
      path_to_seq: path to the sequence which satisfy the assumtion above
      t: ast 
    return
      the updated ast
*)
let var_init_attach : Target.Transfo.local =
  Target.apply_on_path(var_init_attach_aux)


(* [const_non_const_aux t]: This is an auxiliary function for const_non_const
    params:
      t: an ast subterm
    return:
      the updated ast
*)
let const_non_const_aux (t : trm) : trm =
  match t.desc with 
  | Trm_let (vk, (x,tx), init) ->
    begin match vk with
     (* If variable is a constant than whe remove the const and we perform the heap allocation  *)
    | Var_immutable ->  
      trm_let Var_mutable (x, typ_ptr tx) (trm_apps (trm_prim ~loc: t.loc (Prim_new tx)) [init])
    | _ ->
      let var_type = begin match tx.ty_desc with 
      | Typ_ptr t -> t
      | _ -> fail t.loc "const_non_const_aux: expected a pointer type"
      end
      in
      let var_init = begin match init.desc with 
      | Trm_apps(_, [_; init]) -> init
      | _ -> fail t.loc "const_non_const_aux: expected a something of the form 'new ()'"
      end
      in
      trm_let Var_immutable (x,var_type) var_init 
    end
  | _ -> fail t.loc "const_non_const_aux: variable declaration was not matched, make sure the path is correct"

let const_non_const : Target.Transfo.local =
  apply_on_path(const_non_const_aux)


(* [remove_instruction_aux t]: This is an auxiliary function for remove_instruction
    params:
      t: an ast subterm
    return:
      the updated ast
*)
let remove_instruction_aux (_t : trm) : trm =
  (* Replace the current t with an empty sequence *)
  trm_seq ~annot:(Some No_braces) []

let remove_instruction : Target.Transfo.local=
  apply_on_path(remove_instruction_aux)

(* [local_other_name var_type old_var new_var t]: This is an auxiliary function for local_other_name
    params:
      var_type: type of the var
      old_var: old variable for which we want to chang the local name
      new_var: new_variable 
      t: an ast subterm
    return:
      the updated ast
*)
let local_other_name_aux (var_type : typvar) (old_var : var) (new_var : var) (t : trm) : trm =
    begin match t.desc with
    | Trm_seq [no_braces] ->
      begin match no_braces.desc with
        | Trm_seq [f_loop] ->
          begin match f_loop.desc with
          | Trm_for (init, cond, step, body) ->
            let new_type = typ_var var_type (get_typedef var_type) in
            let new_decl = trm_let Var_mutable (new_var, new_type) (trm_apps (trm_prim (Prim_new new_type)) [trm_var old_var])
        
            in
            let new_set_old = trm_set (trm_var old_var) (trm_var new_var) in
            (* let new_del_inst = trm_apps ~annot:(Some Heap_allocated) ~typ:(Some (typ_unit ())) ~is_statement:true (trm_unop (Unop_delete false)) [trm_var new_var] in *)
            let new_loop = trm_for init cond step (change_trm (trm_var old_var)(trm_var new_var) body) in
           
           
              trm_seq (* ~annot:(Some No_braces) *) [
                new_decl;new_loop;new_set_old
              ]
          | _ -> fail t.loc "local_other_name_aux: expected a for loop"
          end
      | _ -> fail t.loc "local_other_name_aux: expected the sequnece which contains the for loop"
      end
    | _ -> fail t.loc "local_other_name_aux: expected the no brace sequence"
    end

let local_other_name (var_type : typvar) (old_var : var) (new_var : var) : Target.Transfo.local = 
  Target.apply_on_path(local_other_name_aux var_type old_var new_var)

(* [delocalize_aux array_size neutral_element fold_operation t]: This is an auxiliary function for deloclize
    params:
      array_size: the size of the array we want to create
      neutral_element: nutral element for reduction phase
      fold_operation: fold_operation for reduction phase
      t: ast subterm
    return:
      the updated ast
*)

(* [insert_trm_aux insert_where t_insert t]: This is an auxiliary function for insert_trm
    params:
      insert_where: a string equal to before or after describing the location where shoudl we insert the trm
      t_insert: the trm to be inserted
      t: the trm around which the trm is going to be inserted
    return:
      the updated ast
*)
let insert_trm_aux (t_insert : trm) (index : int) (t : trm) : trm =
  match t.desc with
  | Trm_seq tl ->
    let tl = Tools.list_insert (index) t_insert tl in
    trm_seq ~annot:t.annot tl
  | _ -> fail t.loc "insert_trm_aux: expected the outer sequence"

let insert_trm (t_insert : trm) (index : int) : Target.Transfo.local =
  Target.apply_on_path(insert_trm_aux t_insert index)