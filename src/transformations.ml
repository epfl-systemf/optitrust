open Ast
open Paths
open Path_constructors
open Translate_ast

let write_log (clog : out_channel) (log : string) : unit =
  output_string clog log; flush clog

(* return the list where the nth element is transformed *)
let change_nth (transfo : 'a -> 'a) (al : 'a list) (n : int) : 'a list =
  List.mapi (fun i a -> if i = n then transfo a else a) al

(* follow an explicit path to apply a function on the corresponding subterm *)
let apply_local_transformation (transfo : trm -> trm) (t : trm)
  (dl : expl_path) : trm =
  let rec aux (dl : expl_path) (t : trm) : trm =
    match dl with
    | [] -> transfo t
    | d :: dl ->
       let annot = t.annot in
       let loc = t.loc in
       let is_instr = t.is_instr in
       let add = t.add in
       let typ = t.typ in
       let attributes = t.attributes in
       begin match d, t.desc with
       | Dir_nth n, Trm_seq tl ->
          trm_seq ~annot ~loc ~add ~attributes (change_nth (aux dl) tl n)
       | Dir_nth n, Trm_array tl ->
          trm_array ~annot ~loc ~add ~typ ~attributes (change_nth (aux dl) tl n)
       | Dir_nth n, Trm_struct tl ->
          trm_struct ~annot ~loc ~add ~typ ~attributes(change_nth (aux dl) tl n)
       | Dir_nth _, Trm_val (Val_array _) ->
          fail loc "apply_local_transformation: val_array should not appear"
       | Dir_nth _, Trm_val (Val_struct _) ->
          fail loc "apply_local_transformation: val_struct should not appear"
       | Dir_cond, Trm_if (cond, then_t, else_t) ->
          trm_if ~annot ~loc ~add ~attributes (aux dl cond) then_t else_t
       | Dir_cond, Trm_while (cond, body) ->
          trm_while ~annot ~loc ~add ~attributes (aux dl cond) body
       | Dir_cond, Trm_for (init, cond, step, body) ->
          trm_for ~annot ~loc ~add ~attributes init (aux dl cond) step body
       | Dir_cond, Trm_switch (cond, cases) ->
          trm_switch ~annot ~loc ~add ~attributes (aux dl cond) cases
       | Dir_then, Trm_if (cond, then_t, else_t) ->
          trm_if ~annot ~loc ~add ~attributes cond (aux dl then_t) else_t
       | Dir_else, Trm_if (cond, then_t, else_t) ->
          trm_if ~annot ~loc ~add ~attributes cond then_t (aux dl else_t)
       | Dir_body, Trm_decl (Def_var (tx, body)) ->
          trm_decl ~annot ~loc ~is_instr ~add ~attributes (Def_var (tx, aux dl body))
       | Dir_body, Trm_decl (Def_fun (x, tx, txl, body)) ->
          trm_decl ~annot ~loc ~is_instr ~add ~attributes (Def_fun (x, tx, txl, aux dl body))
       | Dir_body, Trm_for (init, cond, step, body) ->
          trm_for ~annot ~loc ~add ~attributes init cond step (aux dl body)
       | Dir_body, Trm_while (cond, body) ->
          trm_while ~annot ~loc ~add ~attributes cond (aux dl body)
       | Dir_body, Trm_abort (Ret (Some body)) ->
          trm_abort ~annot ~loc ~add ~attributes (Ret (Some (aux dl body)))
       | Dir_body, Trm_labelled (l, body) ->
          trm_labelled ~annot ~loc ~add ~attributes l (aux dl body)
       | Dir_body, Trm_decoration(left, body, right) ->
          trm_decoration ~annot ~loc ~add ~attributes left right (aux dl body) 
       | Dir_for_init, Trm_for (init, cond, step, body) ->
          trm_for ~annot ~loc ~add ~attributes (aux dl init) cond step body
       | Dir_for_step, Trm_for (init, cond, step, body) ->
          trm_for ~annot ~loc ~add ~attributes init cond (aux dl step) body
       | Dir_app_fun, Trm_apps (f, tl) ->
          (*
            warning: the type of f may change
            -> print and reparse to have the right type
           *)
          trm_apps ~annot ~loc ~is_instr ~add ~typ ~attributes (aux dl f) tl
       | Dir_arg n, Trm_apps (f, tl) ->
          trm_apps ~annot ~loc ~is_instr ~add ~typ ~attributes f
            (change_nth (aux dl) tl n)
       | Dir_arg n, Trm_decl (Def_fun (x, tx, txl, body)) ->
          let txl' =
            change_nth
              (fun (x, tx) ->
                let t' = aux dl (trm_var ~loc x) in
                match t'.desc with
                | Trm_var x' -> (x', tx)
                | _ ->
                   fail loc ("apply_local_transformation: transformation " ^
                               "must preserve fun arguments")
              )
              txl
              n
          in
          trm_decl ~annot ~loc ~is_instr ~add ~attributes
            (Def_fun (x, tx, txl', body))
       | Dir_name, Trm_decl (Def_var ((x, tx), body)) ->
          let t' = aux dl (trm_var ~loc x) in
          begin match t'.desc with
          | Trm_var x' ->
             trm_decl ~annot ~loc ~is_instr ~add ~attributes
               (Def_var ((x', tx), body))
          | _ ->
             fail loc ("apply_local_transformation: transformation " ^
                         "must preserve names")
          end
       | Dir_name, Trm_decl (Def_fun (x, tx, txl, body)) ->
          let t' = aux dl (trm_var ~loc x) in
          begin match t'.desc with
          | Trm_var x' ->
             trm_decl ~annot ~loc ~is_instr ~add ~attributes
               (Def_fun (x', tx, txl, body))
          | _ ->
             fail loc ("apply_local_transformation: transformation " ^
                         "must preserve names")
          end
       | Dir_name, Trm_labelled (l, body) ->
          let t' = aux dl (trm_var ~loc l) in
          begin match t'.desc with
          | Trm_var l' ->
             trm_labelled ~annot ~loc ~add ~attributes l' body
          | _ ->
             fail loc ("apply_local_transformation: transformation " ^
                         "must preserve names")
          end
          
       | Dir_case (n, cd), Trm_switch (cond, cases) ->
          trm_switch ~annot ~loc ~add ~attributes cond
            (change_nth
               (fun (tl, body) ->
                 match cd with
                 | Case_body -> (tl, aux dl body)
                 | Case_name i ->
                    (change_nth (fun ith_t -> aux dl ith_t) tl i, body)
               )
               cases
               n
            )
        | _, _ ->
           let s = string_of_dir d in
           fail loc ("apply_local_transformation: direction " ^ s ^
                       " does not match")
       end

  in
  aux dl t

(* insert a after rank n in the list *)
let rec list_insert (n : int) (a : 'a) (al : 'a list) : 'a list =
  if n < 0 then a :: al else List.hd al :: list_insert (n - 1) a (List.tl al)

(*
  insert inert after the subterm pointed at by dl in t
  assumption: dl points at a seq element, thus ends with Dir_nth n
  if the inserted element must be first in the seq, use n < 0
 *)
let insert_trm_after (dl : expl_path) (insert : trm) (t : trm) : trm =
  let dl' = List.rev dl in
  match List.hd dl' with
  | Dir_nth n ->
     apply_local_transformation
       (fun t' ->
         match t'.desc with
         | Trm_seq tl ->
            trm_seq ~annot:t'.annot ~loc:t'.loc ~add:t'.add
              ~attributes:t'.attributes (list_insert n insert tl)
         | _ -> fail t'.loc "insert_trm_after: path points at wrong term"
       )
       t
       (List.rev (List.tl dl'))
  | _ -> fail t.loc "insert_trm_after: bad path"


let add_label (label : string) (pl : path list) (t : trm) : trm =
  let p = List.flatten pl in
  let b = !Flags.verbose in
  Flags.verbose := false;
  let epl = resolve_path p t in
  Flags.verbose := b;
  match epl with
  | [] ->
     print_info t.loc "add_label: no matching subterm\n";
     t
  | [dl] -> apply_local_transformation (trm_labelled label) t dl
  | _ ->
     (*
         folding works since no path in epl is the prefix of a subsequent path
      *)
     foldi
       (fun i -> apply_local_transformation
                   (trm_labelled (label ^ "_" ^ string_of_int i)))
       t epl



let list_remove x xs = List.filter (fun y -> y <> x) xs 


let list_remove_set ys xs = List.fold_left (fun acc y -> list_remove y acc) xs ys 

let move_fields_before x local_l l = 
let l = list_remove_set l local_l in 
let rec aux acc = function 
| [] -> acc (* raise an error x not part of the list *)
| hd :: tl -> if hd = x then aux (local_l @ hd :: acc) tl (* local_l @ hd :: acc @ tl *)
else aux (hd :: acc) tl 
in aux [] l

(* 
  - tail recursive approach => more efficient 
  - non-tail rec => easier to read 

     let rec insert_after x xs l =
        match l with
        | [] -> error
        | y::q -> if x = y then xs@l else y::(insert_after x xs q)
*)




let move_fields_after x local_l l = 
  let l = list_remove_set l local_l in 
  let rec aux acc = function 
    | [] -> acc
    | hd :: tl -> 
        if hd = x
          then aux (hd :: local_l @ acc) tl 
          else aux (hd :: acc) tl 
    in
  aux [] l


let fields_reorder_aux (clog :out_channel) ?(struct_fields : fields = []) ?(move_before : field = "") ?(move_after : field = "")(t : trm) : trm  = 
    let log : string = 
      let loc : string = 
        match t.loc with 
        | None -> ""
        | Some (_, line) -> Printf.sprintf "at line %d " line
      in Printf.sprintf
          ("  - expression\n%s\n" ^^
          "    %sis a struct type\n"
          )
      (ast_to_string t) loc 
    in
    write_log clog log;
    begin match t.desc with
      | Trm_decl (Def_typ (x,dx)) ->
        
        let field_list, field_map = 
          match dx.ty_desc with
            | Typ_struct(l,m,_) -> l,m
            |_ -> fail t.loc "fields_reorder: the type should be a typedef struct"
          in
        let reordered_fields = 
          match move_before, move_after with 
          | "",_ -> move_fields_after move_after struct_fields field_list
          | _, "" -> move_fields_before move_before struct_fields field_list
          | _,_-> fail t.loc "fields_reorder: only one of move_before or move_after should be specified"
          in
          (* TODO: field order should match in the AST what works. *)
        let t_yp = {ty_desc = Typ_struct(List.rev reordered_fields,field_map,x); ty_annot = dx.ty_annot; ty_attributes = dx.ty_attributes} in
        
        trm_decl ~annot:t.annot ~loc:t.loc ~is_instr:t.is_instr ~add:t.add
          ~attributes:t.attributes (Def_typ(x,t_yp) )    

      | _ -> fail t.loc "fields_reorder: expected a definiton"
      end
    
 

let fields_reorder (clog :out_channel) ?(struct_fields : fields = []) ?(move_before : field = "") ?(move_after : field = "") (pl : path list) (t : trm) : trm  = 
  let p = List.flatten pl in 
  let b = !Flags.verbose in
  Flags.verbose := false;
  let epl = resolve_path p t in 
  Flags.verbose := b;
  match epl with 
  | [] -> 
      print_info t.loc "Struct field reordering\n";
      t
  | _ -> 
      List.fold_left 
        (fun t dl -> 
          apply_local_transformation (fields_reorder_aux clog ~struct_fields ~move_before ~move_after) t dl )
        t
        epl


 

 
 

let left_decoration (index:int):string  = "/*@" ^ string_of_int index ^ "<*/"  

let right_decoration (index:int):string  = "/*>" ^ string_of_int index ^ "@*/"
   

let  show_path(pl : path list)  (t : trm) : trm =
  let p = List.flatten pl in 
  let b = !Flags.verbose in 
  Flags.verbose := false;
  let epl = resolve_path p t in 
  Flags.verbose := b;
  match epl with 
  | [] ->
    print_info t.loc "show_path: not matching subterm\n";
    t
  | [dl] -> apply_local_transformation (trm_decoration (left_decoration 0) (right_decoration 0)) t dl
  
  | _ ->
     (*
         folding works since no path in epl is the prefix of a subsequent path
      *)
     foldi
       (fun i -> apply_local_transformation
                   (trm_decoration (left_decoration i) (right_decoration i )))
       t epl


let rec delete_label (label : string) (t : trm) : trm =
  match t.desc with
  | Trm_labelled (l, t') when l = label -> t'
  | _ -> trm_map (delete_label label) t

let rec delete_path_decorators (t : trm) : trm = 
  match t.desc with 
  | Trm_decoration (_,t',_) -> t'
  | _ -> trm_map (delete_path_decorators ) t




(* delete the labels which have a prefix in the list *)
let delete_labels (sl : string list) (t : trm) : trm =
  let rec aux (s : string) (t : trm) : trm =
    match t.desc with
    | Trm_labelled (l, t')
         when Str.string_match (Str.regexp (Str.quote s)) l 0 ->
       t'
    | _ -> trm_map (aux s) t
  in
  List.fold_left (fun t l -> aux l t) t sl

(* make sure each occurence of y in t is marked with type variable x *)
let rec replace_type_with (x : typvar) (y : var) (t : trm) : trm =
  match t.desc with
  | Trm_var y' when y' = y ->
     trm_var ~annot:t.annot ~loc:t.loc ~add:t.add ~typ:(Some (typ_var x))
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
  let rec aux (t : trm) : ilset funmap =
    match t.desc with
    | Trm_decl (Def_var (_, body)) -> aux body
    | Trm_decl (Def_fun (_, _, _, body)) -> aux body
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
                | Some {ty_desc = Typ_var x'; _} when x' = x -> i :: il
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
          match path_to_decl f global_trm with
          (* if the declaration cannot be found, ignore this function *)
          | None ->
             print_info global_trm.loc
               ("functions_with_arg_type: cannot find declaration of " ^^
                  "function %s, ignoring it.\n") f;
             Fun_map.remove f res
          | Some dl ->
             let (def, _) = resolve_explicit_path dl global_trm in
             begin match def.desc with
             | Trm_decl (Def_fun (_, _, args, body)) ->
                let b = replace_arg_types_with x il args body in
                (* then do a recursive call on the new body *)
                res +@ functions_with_arg_type ~outer_trm:(Some global_trm) x b
             | _ ->
                fail t.loc
                  ("functions_with_arg_type: wrong path to declaration of " ^ f)
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
      condition on annotation: toplevel declarations might contain a heap
      allocated variable and hence we can find a no_brace seq inside a
      delete_instructions seq, which we do not want to inline
     *)
    | Trm_seq tl when t.annot <> Some Delete_instructions ->
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
         match path_to_decl f t' with
         | None ->
            fail t'.loc
              ("insert_fun_copies: cannot find declaration of function " ^ f)
         | Some dl ->
            let (fdecl, _) = resolve_explicit_path dl t' in
            begin match fdecl.desc with
            | Trm_decl (Def_fun (f', r, tvl, b)) when f = f' ->
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
                         (change_nth (fun (y, _) -> (y, typ_var x))) tvl il
                     in
                     (* add index to labels in the body of the function *)
                     let b' =
                       label_aux i (replace_arg_types_with x il tvl' b)
                     in
                     (* create the copy of f corresponding to il *)
                     (trm_decl
                        (Def_fun (name f ^ "_" ^ string_of_int i, r, tvl', b'))
                     ) ::
                     tl
                   )
                   ils
                   []
               in
               (* insert the copies of f *)
               insert_trm_after dl
                 (trm_seq ~annot:(Some No_braces) (List.rev tl)) t'
            | _ -> fail t'.loc "insert_fun_copies: bad path to fun decl"
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
  let is_instr = t.is_instr in
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
                    | Some {ty_desc = Typ_var x'; _} when x' = x -> i :: il
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
              trm_apps ~annot ~loc ~is_instr ~add ~typ ~attributes
                (trm_var ~annot:fun_.annot ~loc:fun_.loc ~add:fun_.add
                   ~attributes:fun_.attributes f') args'
           end
        end
     | _ -> trm_map (replace_fun_names name ilsm x) t
     end
  | _ -> trm_map (replace_fun_names name ilsm x) t

(*
  insert t_inserted either before the position pointed at by insert_before or
  after the position pointed at by insert_after in t
  both must be resolved as paths to a seq element
 *)
let insert_trm ?(insert_before : path list = [])
  ?(insert_after : path list = []) (t_inserted : trm) (t : trm) : trm =
  let p =
    match insert_before, insert_after with
    | [], _ :: _ -> List.flatten insert_after
    | _ :: _, [] -> List.flatten insert_before
    | [], [] -> fail t.loc "insert_trm: please specify an insertion point"
    | _ -> fail t.loc "insert_trm: cannot insert both before and after"
  in
  let b = !Flags.verbose in
  Flags.verbose := false;
  let epl = resolve_path p t in
  Flags.verbose := b;
  match epl with
  | [] ->
     print_info t.loc "insert_trm: no matching subterm\n";
     t
  | _ ->
     List.fold_left
       (fun t' dl ->
         match List.rev dl with
         | Dir_nth n :: dl' ->
            begin match insert_before, insert_after with
            (* insert after *)
            | [], _ :: _ -> insert_trm_after dl t_inserted t'
            (* insert before: replace n with n - 1 *)
            | _ :: _, [] ->
               insert_trm_after (List.rev (Dir_nth (n - 1) :: dl')) t_inserted
                 t'
            | [], [] ->
               fail t'.loc "insert_trm: please specify an insertion point"
            | _ -> fail t'.loc "insert_trm: cannot insert both before and after"
            end
         | _ -> fail t'.loc "insert_trm: bad insertion path"
       )
       t
       epl

(*
  replace occurrences of t_before with t_after in t
  paths point at subterms in which all occurences will be replaced
  the empty path means all occurences will be replaced (default behaviour)
  assumption: t_before and t_after are equivalent (in terms of value and of side
  effects)
 *)
let change_trm ?(change_at : path list list = [[]]) (t_before : trm)
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
           when t'.annot = Some Heap_allocated ->
         trm_seq ~annot:t'.annot ~loc:t'.loc ~add:t'.add
           ~attributes:t'.attributes
           [
             t_decl;
             trm_set ~annot:(Some Initialisation_instruction) ~loc lhs
               (apply_change init)
           ]
      | _ -> trm_map apply_change t'
  in
  List.fold_left
    (fun t' pl ->
      let p = List.flatten pl in
      let b = !Flags.verbose in
      Flags.verbose := false;
      let epl = resolve_path p t' in
      Flags.verbose := b;
      match epl with
      | [] ->
         print_info t'.loc "change_trm: no matching subterm for path %s\n"
           (string_of_path p);
         t'
      | _ -> List.fold_left (apply_local_transformation apply_change) t' epl
    )
    t
    change_at


(* same as change_trm but for types *)
let change_typ ?(change_at : path list list = [[]]) (ty_before : typ)
  (ty_after : typ) (t : trm) : trm =
  (* change all occurences of ty_before in ty *)
  let rec change_typ (ty : typ) : typ =
    (* necessary because of annotations in trms that may be different *)
    if typ_to_string ty = typ_to_string ty_before then ty_after
    else typ_map change_typ ty
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
      | Trm_decl (Def_var ((y, ty), init)) ->
         trm_decl ~annot:t.annot ~loc:t.loc ~is_instr:t.is_instr ~add:t.add
           ~attributes:t.attributes (Def_var ((y, change_typ ty), aux init))
      | Trm_decl (Def_fun (f, ty, args, body)) ->
         trm_decl ~annot:t.annot ~loc:t.loc ~is_instr:t.is_instr ~add:t.add
           ~attributes:t.attributes
           (Def_fun (f, change_typ ty,
                     List.map (fun (y, ty) -> (y, change_typ ty)) args,
                     aux body)
           )
      | Trm_decl (Def_typ (y, ty)) ->
         trm_decl ~annot:t.annot ~loc:t.loc ~is_instr:t.is_instr ~add:t.add
            ~attributes:t.attributes (Def_typ (y, change_typ ty))
      | _ -> trm_map aux t
    in
    replace_type_annot (aux t)
  in
  List.fold_left
    (fun t' pl ->
      let p = List.flatten pl in
      let b = !Flags.verbose in
      Flags.verbose := false;
      let epl = resolve_path p t' in
      Flags.verbose := b;
      match epl with
      | [] ->
         print_info t'.loc "change_typ: no matching subterm for path %s\n"
           (string_of_path p);
         t'
      | _ -> List.fold_left (apply_local_transformation apply_change) t' epl
    )
    t
  
    change_at



(*
  find the definition x = dx pointed at by pl and replace occurrences of dx with
  x
  paths point at subterms in which all occurences will be replaced
  the empty path means all occurences will be replaced (default behaviour)
  as_reference option for variable declarations: if dx = &dx' replace dx' with
  *x instead of &dx' with x
 *)
let fold_decl (clog : out_channel) ?(as_reference : bool = false)
  ?(fold_at : path list list = [[]]) (pl : path list) (t : trm) : trm =
  let p = List.flatten pl in
  let b = !Flags.verbose in
  Flags.verbose := false;
  let epl = resolve_path p t in
  Flags.verbose := b;
  match epl with
  | [dl] ->
     let (t_def, _) = resolve_explicit_path dl t in
     let log : string =
       Printf.sprintf
         ("  - expression\n%s\n" ^^
          if as_reference then
          "    is a variable declaration of the form\n" ^^
          "      type* x = &dx\n"
          else
          "    is a variable/type declaration\n"
         )
         (ast_to_string t_def)
     in
     write_log clog log;
     begin match t_def.desc with
     (* const variables *)
     | Trm_decl (Def_var ((x, _), dx)) ->
        let t_x =
          if as_reference then trm_apps (trm_unop Unop_get) [trm_var x]
          else trm_var x
        in
        let def_x =
          if not as_reference then dx
          else
            match dx.add with
            | Add_address_of_operator :: addl -> {dx with add = addl}
            | _ -> fail t_def.loc "fold_decl: expected a reference"
        in
        let t = change_trm ~change_at:fold_at def_x t_x t in
        (*
          def_x might have been replaced with x in the definition of x
          -> replace it again with def_x
         *)
        let change_at =
          [[cVarDef ~name:x ~body:[cVar ~name:x ()] (); cBody ~strict:true ()]]
        in
        change_trm ~change_at t_x def_x t
     (*
       heap allocated variables
       note: an initialisation must be given
      *)
     | Trm_seq [{desc = Trm_decl (Def_var ((x, _), _)); _};
                {desc = Trm_apps (_, [_; dx]); _}]
          when t_def.annot = Some Heap_allocated ->
        let t_x =
          trm_apps ~annot:(Some Heap_allocated) (trm_unop Unop_get) [trm_var x]
        in
        let t_x =
          if as_reference then trm_apps (trm_unop Unop_get) [t_x] else t_x
        in
        let def_x =
          if not as_reference then dx
          else
            match dx.add with
            | Add_address_of_operator :: addl -> {dx with add = addl}
            | _ -> fail t_def.loc "fold_decl: expected a reference"
        in
        let t = change_trm ~change_at:fold_at def_x t_x t in
        (* make sure def_x is not replaced in the definition of x here too *)
        let change_at =
          [[cVarDef ~name:x ~body:[cVar ~name:x ()] (); cNth ~strict:true 1;
            cArg ~strict:true 1]]
        in
        change_trm ~change_at t_x def_x t
     (* typedef *)
     | Trm_decl (Def_typ (x, dx)) ->
        let ty_x = typ_var x in
        let t = change_typ ~change_at:fold_at dx ty_x t in
        (* make sure dx is not replaced in the definition of x here too *)
        let change_at = [[cType ~name:x ()]] in
        change_typ ~change_at ty_x dx t
     (* fun decl *)
     | Trm_decl (Def_fun _) ->
        fail t.loc "fold_decl: fun declaration folding is unsupported"
     | _ -> fail t.loc "fold_decl: expected a definition"
     end
  | _ -> fail t.loc "fold_decl: the path must point to exactly 1 subterm"

(*
  insert a definition x = dx either before the position pointed at by
  insert_before or after the position pointed at by insert_after
  both must be resolved as paths to a seq element
  x may be a const variable or not (not const by default)
  option: make x a reference (x = &dx)
  assumptions:
    - no conflicts with the new name x
    - for a given seq, the insertion path points to at most one of its elements
    - if x is a reference, dx denotes a memory cell
 *)
let insert_decl ?(insert_before : path list = [])
  ?(insert_after : path list = []) ?(const : bool = false)
  ?(as_reference : bool = false) (x : var) (dx : trm) (t : trm) : trm =
  let tx =
    match dx.typ with
    | None -> fail dx.loc "insert_decl: cannot find definition type"
    | Some tx -> if as_reference then typ_ptr tx else tx
  in
  let def_x =
    if as_reference then {dx with add = Add_address_of_operator :: dx.add}
    else dx
  in
  let t_insert =
    if const then trm_decl (Def_var ((x, tx), def_x))
    else
      trm_seq ~annot:(Some Heap_allocated)
        [trm_decl (Def_var ((x, typ_ptr tx), trm_prim (Prim_new tx)));
         trm_set ~annot:(Some Initialisation_instruction) (trm_var x) def_x
        ]
  in
  (* compute the explicit path for later use *)
  let p =
    match insert_before, insert_after with
    | [], _ :: _ -> List.flatten insert_after
    | _ :: _, [] -> List.flatten insert_before
    | [], [] -> fail t.loc "insert_decl: please specify an insertion point"
    | _ -> fail t.loc "insert_decl: cannot insert both before and after"
  in
  let b = !Flags.verbose in
  Flags.verbose := false;
  let epl = resolve_path p t in
  Flags.verbose := b;
  (* insert the definition *)
  let t = insert_trm ~insert_before ~insert_after t_insert t in
  (*
    don't forget the delete instruction if x is heap allocated
    use explicit path because p will not be resolved as the position of the
    definition
   *)
  if const then t
  else
    (*
      add a seq with delete instruction around the pointed term containing the
      declaration
     *)
    let create_delete_instr (dl : expl_path) (t : trm) : trm =
      apply_local_transformation
        (fun t ->
          (* t is expected to be a seq *)
          trm_seq ~annot:(Some Delete_instructions)
            [t;
             trm_apps ~annot:(Some Heap_allocated)
               ~typ:(Some (typ_unit ()))
               (trm_unop (Unop_delete false)) [trm_var x]
            ]
        )
        t
        dl
    in
    List.fold_left
      (fun t dl ->
        match List.rev dl with
        (*
          the seq containing the definition might be inside a seq with delete
          instructions
          -> do not create a seq
         *)
        | Dir_nth _ :: Dir_nth n :: dl ->
           apply_local_transformation
             (fun t ->
               match t.desc with
               | Trm_seq (t' :: del_instr_l)
                    when t.annot = Some Delete_instructions ->
                  trm_seq ~annot:(Some Delete_instructions)
                    (t' ::
                     (trm_apps ~annot:(Some Heap_allocated)
                        ~typ:(Some (typ_unit ()))
                        (trm_unop (Unop_delete false)) [trm_var x]) ::
                     del_instr_l
                    )
               (*
                 if we do not find a seq of delete instructions, go deeper to
                 create the seq
                *)
               | _ -> create_delete_instr [Dir_nth n] t
             )
             t
             (List.rev dl)
        | Dir_nth _ :: dl -> create_delete_instr (List.rev dl) t
        | _ -> fail t.loc "insert_definition: expected a path to a seq"
      )
      t
      epl

(* same as insert_definition but for a constant *)
let insert_const ?(insert_before : path list = [])
  ?(insert_after : path list = []) (x : var) (dx : trm) (t : trm) : trm =
  insert_decl ~insert_before ~insert_after ~const:true x dx t

(*
  insert a type declaration x = dx either before the position pointed at by
  insert_before or after the position pointed at by insert_after
  both must be resolved as paths to a seq element
  assumption: no conflicts with the new name x
 *)
let insert_typedef ?(insert_before : path list = [])
  ?(insert_after : path list = []) (x : typvar) (dx : typ) (t : trm) : trm =
  insert_trm ~insert_before ~insert_after (trm_decl (Def_typ (x, dx))) t

(*
  combine insert_definition and fold_decl
  assumption: if x is not a reference, no effects for dx and it has the same
  value through all its occurences
 *)
let insert_and_fold (clog : out_channel) ?(insert_before : path list = [])
  ?(insert_after : path list = []) ?(const : bool = false)
  ?(as_reference : bool = false) ?(fold_at : path list list = [[]]) (x : var)
  (dx : trm) (t : trm) : trm =
  (* compute the explicit path for later use *)
  let p =
    match insert_before, insert_after with
    | [], _ :: _ -> List.flatten insert_after
    | _ :: _, [] -> List.flatten insert_before
    | [], [] -> fail t.loc "insert_and_fold: please specify an insertion point"
    | _ -> fail t.loc "insert_and_fold: cannot insert both before and after"
  in
  let b = !Flags.verbose in
  Flags.verbose := false;
  let epl = resolve_path p t in
  Flags.verbose := b;
  (* insert the definition *)
  let t =
    insert_decl ~insert_before ~insert_after ~const ~as_reference x dx t
  in
  (*
    fold the definition
    use explicit path because p will not be resolved as the position of the
    definition
    any path in epl is ok to do so
   *)
  match epl with
  | [] -> fail t.loc "insert_and_fold: no insertion point"
  | dl :: _ ->
     let def_pathl =
       let pathl_of_expl_path (dl : expl_path) : path list =
         List.map (fun d -> [Constr_strict; Constr_dir d]) dl
       in
       match List.rev dl with
       | Dir_nth n :: dl ->
          let n =
            match insert_before, insert_after with
            (* insert after: add 1 to n *)
            | [], _ :: _ -> n + 1
            (* insert before: n is the position of the definition *)
            | _ :: _, [] -> n
            | [], [] ->
               fail t.loc "insert_and_fold: please specify an insertion point"
            | _ ->
               fail t.loc "insert_and_fold: cannot insert both before and after"
          in
          let (t_container, _) = resolve_explicit_path (List.rev dl) t in
          begin match t_container.annot with
          (*
            in case of heap allocation, a seq (for delete instructions) may be
            added around the last container
            -> add a nth 0 direction before the last direction if it is the case
           *)
          | Some Delete_instructions ->
             pathl_of_expl_path (List.rev (Dir_nth n :: Dir_nth 0 :: dl))
          | _ -> pathl_of_expl_path (List.rev (Dir_nth n :: dl))
          end
       | _ -> fail t.loc "insert_and_fold: expected a path to a seq element"
     in
     (* replace dx with &dx before folding if we have a reference *)
     fold_decl clog ~as_reference ~fold_at def_pathl t

(* same as insert_and_fold but for types *)
let insert_and_fold_typedef (clog : out_channel)
  ?(insert_before : path list = []) ?(insert_after : path list = [])
  ?(fold_at : path list list = [[]]) (x : typvar) (dx : typ) (t : trm) : trm =
  (* compute the explicit path for later use *)
  let p =
    match insert_before, insert_after with
    | [], _ :: _ -> List.flatten insert_after
    | _ :: _, [] -> List.flatten insert_before
    | [], [] ->
       fail t.loc "insert_and_fold_typedef: please specify an insertion point"
    | _ ->
       fail t.loc "insert_and_fold_typedef: cannot insert both before and after"
  in
  let b = !Flags.verbose in
  Flags.verbose := false;
  let epl = resolve_path p t in
  Flags.verbose := b;
  (* insert the typedef *)
  let t = insert_typedef ~insert_before ~insert_after x dx t in
  (*
    fold the typedef
    use explicit path because p will not be resolved as the position of the
    definition
    any path in epl is ok to do so
   *)
  match epl with
  | [] -> fail t.loc "insert_and_fold_typedef: no insertion point"
  | dl :: _ ->
     let dl =
       match List.rev dl with
       | Dir_nth n :: dl' ->
          let n =
            match insert_before, insert_after with
            (* insert after: add 1 to n *)
            | [], _ :: _ -> n + 1
            (* insert before: n is the position of the definition *)
            | _ :: _, [] -> n
            | [], [] ->
               fail t.loc
                 "insert_and_fold_typedef: please specify an insertion point"
            | _ ->
               fail t.loc
                 "insert_and_fold_typedef: cannot insert both before and after"
          in
          List.rev (Dir_nth n :: dl')
       | _ -> fail t.loc "insert_and_fold_typedef: expected a path to a seq"
     in
     let def_pathl = List.map (fun d -> [Constr_strict; Constr_dir d]) dl in
     fold_decl clog ~fold_at def_pathl t

let filteri (f : int -> 'a -> bool) (al : 'a list) : 'a list =
  let aol = List.mapi (fun i a -> if f i a then Some a else None) al in
  List.filter_map (fun ao -> ao) aol

(*
  remove the declaration pointed at by pl
  pl must be resolved as a path to a seq element
  assumption: the declared object is not used in t
 *)
let remove_decl (clog : out_channel) (pl : path list) (t : trm) : trm =
  let p = List.flatten pl in
  let b = !Flags.verbose in
  Flags.verbose := false;
  let epl = resolve_path p t in
  Flags.verbose := b;
  match epl with
  | [dl] ->
     (* get the declaration for later use *)
     let (t_decl, _) = resolve_explicit_path dl t in
     let log : string =
       let loc : string =
         match t_decl.loc with
         | None -> ""
         | Some (_, line) -> Printf.sprintf "at line %d " line
       in
       Printf.sprintf
         ("  - expression\n%s\n" ^^
          "    %sis a declaration\n"
         )
         (ast_to_string t_decl) loc
     in
     write_log clog log;
     let log : string =
       let x = decl_name t_decl in
       Printf.sprintf "  - %s is not used in the remainder of the program\n" x
     in
     write_log clog log;
     let dl = List.rev dl in
     let n =
       match List.nth_opt dl 0 with
       | Some (Dir_nth n) -> n
       | _ -> fail t.loc "remove_decl: the path must point at a seq element"
     in
     let t =
       apply_local_transformation
         (fun (t : trm) ->
           match t.desc with
           | Trm_seq tl ->
              let tl = filteri (fun i _ -> i <> n) tl in
              trm_seq ~annot:t.annot ~loc:t.loc ~add:t.add tl
           | _ -> fail t.loc "remove_decl: expected a seq"
         )
         t
         (* remove the last direction to point at the seq *)
         (List.rev (List.tl dl))
     in
     (* remove delete instruction if the declaration is a heap allocation *)
     begin match t_decl.desc with
     | Trm_seq _ when t_decl.annot = Some Heap_allocated ->
        let x = decl_name t_decl in
        apply_local_transformation
          (fun (t : trm) ->
            match t.desc with
            | Trm_seq (t_body :: del_instr_l)
                 when t.annot = Some Delete_instructions ->
               let del_instr_l =
                 List.filter
                   (fun (t_del : trm) ->
                     match t_del.desc with
                     | Trm_apps (_, [{desc = Trm_var y; _}]) when y = x -> false
                     | _ -> true
                   )
                   del_instr_l
               in
               begin match del_instr_l with
               | [] -> t_body
               | _ ->
                  trm_seq ~annot:t.annot ~loc:t.loc ~add:t.add
                    ~attributes:t.attributes (t_body :: del_instr_l)
               end
            | _ -> fail t.loc "remove_decl: expected delete instructions"
          )
          t
          (*
            remove the two last directions to point at the seq containing the
            delete instructions
           *)
          (List.rev (List.tl (List.tl dl)))
     | _ -> t
     end
  | _ -> fail t.loc "remove_decl: the path must point at exactly 1 subterm"

(* compute a fresh variable (w.r.t. t) based on x *)
let fresh_in (t : trm) (x : var) : var =
  if not (is_used_var_in t x) then x
  else
    begin
      let n = ref 0 in
      while is_used_var_in t (x ^ "_" ^ string_of_int !n) do
        incr n
      done;
      x ^ "_" ^ string_of_int !n
    end

let eliminate_goto_next (t : trm) : trm =
  let rec elim_in_list (tl : trm list) : trm list =
    match tl with
    | t1 :: t2 :: tl ->
       begin match t1.desc, t2.desc with
       | Trm_goto l1, Trm_labelled (l2, _) when l1 = l2 ->
          elim_in_list (t2 :: tl)
       | _ -> t1 :: (elim_in_list (t2 :: tl))
       end
    | _ -> tl
  in
  let rec aux (t : trm) : trm =
    match t.desc with
    | Trm_seq tl ->
       trm_seq ~annot:t.annot ~loc:t.loc ~add:t.add ~attributes:t.attributes
         (elim_in_list (List.map aux tl))
    | _ -> trm_map aux t
  in
  aux t

let group_decl_init (t : trm) : trm =
  let rec group_in_list (tl : trm list) : trm list =
    match tl with
    | t1 :: t2 :: tl ->
       begin match t1.desc, t2.desc with
       | Trm_seq [{desc = Trm_decl (Def_var ((x, tx), dx)); _}],
         Trm_apps ({desc = Trm_val (Val_prim (Prim_binop Binop_set)); _},
                   [{desc = Trm_var y; _}; init])
             when y = x && t1.annot = Some Heap_allocated ->
          let t =
            trm_seq ~annot:(Some Heap_allocated) ~loc:t1.loc
              [
                trm_decl (Def_var ((x, tx), dx));
                trm_set ~annot:(Some Initialisation_instruction) (trm_var x)
                  init
              ]
          in
          group_in_list (t :: tl)
       | _ -> t1 :: (group_in_list (t2 :: tl))
       end
    | _ -> tl
  in
  let rec aux (t : trm) : trm =
    match t.desc with
    | Trm_seq tl ->
       trm_seq ~annot:t.annot ~loc:t.loc ~add:t.add ~attributes:t.attributes
         (group_in_list (List.map aux tl))
    | _ -> trm_map aux t
  in
  aux t

(*
  instr containing f(arg1, …, argn) is replaced with
  {
    x1 = arg1
    …
    xn = argn
    decl result
    body[x1, …, xn][return r := {result = r; goto return_label}]
    return_label:
      instr[f(arg1, …, argn) := result]
  }
  if tf is void, result won't be used, but instead the empty statement
 *)
let inline_fun_decl ?(inline_at : path list list = [[]]) (result : var)
  (return_label : label) (f : var) (tf : typ) (args : typed_var list)
  (body : trm) (t : trm) : trm =
  (* new names replacing the argument names *)
  let fresh_args = List.map (fun (x, tx) -> (fresh_in t x, tx)) args in
  (* name for the result of f, result might be an argument name *)
  let result =
    fresh_in
      (trm_seq
         (t ::
            List.map
              (fun x_tx -> trm_decl (Def_var (x_tx, trm_lit Lit_uninitialized)))
              fresh_args
         )
      )
      result
  in
  (* result is heap allocated *)
  let result_decl =
    trm_seq ~annot:(Some Heap_allocated)
      [trm_decl (Def_var ((result, typ_ptr tf), trm_prim (Prim_new tf)))]
  in
  (* body where the argument names are substituted *)
  let body =
    List.fold_left
      (fun body (x, _) ->
        change_trm
          (trm_var x)
          (* arguments will be heap allocated *)
          (trm_apps ~annot:(Some Heap_allocated) (trm_unop Unop_get)
             [trm_var (fresh_in t x)])
          body
      )
      body
      args
  in
  (* body where res is used instead of return statements *)
  let replace_return (t : trm) : trm =
    let rec aux (t : trm) : trm =
      match t.desc with
      (* remove delete instruction related to return statement if any *)
      | Trm_seq tl when t.annot = Some Delete_instructions ->
         begin match List.rev tl with
         | {desc = Trm_abort (Ret (Some r)); _} :: _ ->
            trm_seq ~annot:(Some No_braces) ~loc:t.loc
              [trm_set ~loc:t.loc ~is_instr:true (trm_var result) r;
               trm_goto ~loc:t.loc return_label]
         | {desc = Trm_abort (Ret None); _} :: _ ->
            trm_goto ~loc:t.loc return_label
         | _ -> trm_map aux t
         end
      | Trm_abort (Ret (Some r)) ->
         trm_seq ~annot:(Some No_braces) ~loc:t.loc
           [trm_set ~loc:t.loc ~is_instr:true (trm_var result) r;
            trm_goto ~loc:t.loc return_label]
      | Trm_abort (Ret None) -> trm_goto ~loc:t.loc return_label
      | _ -> trm_map aux t
    in
    clean_up_no_brace_seq (aux t)
  in
  let body = replace_return body in
  let bodyl =
    match body.annot with
    | Some Delete_instructions -> [body]
    | _ ->
       begin match body.desc with
       | Trm_seq tl -> tl
       | _ -> [body]
       end
  in
  (* inline f everywhere in t *)
  let rec apply_change (t : trm) : trm =
    (* we look for instructions that contain a call to f *)
    if not (t.is_instr && contains_call_to_fun f t)
    then trm_map apply_change t
    else
      let arg_vals = fun_call_args f t in
      let arg_decls =
        List.map2
          (fun (x, tx) dx ->
            trm_seq ~annot:(Some Heap_allocated)
              [
                trm_decl (Def_var ((x, typ_ptr tx), trm_prim (Prim_new tx)));
                trm_set ~annot:(Some Initialisation_instruction) (trm_var x) dx
              ]
          )
          fresh_args
          arg_vals
      in
      let arg_dels =
        List.rev_map
          (fun (x, _) ->
            trm_apps ~annot:(Some Heap_allocated) ~typ:(Some (typ_unit ()))
              (trm_unop (Unop_delete false)) [trm_var x]
          )
          fresh_args
      in
      let t =
        match tf.ty_desc with
        | Typ_unit ->
           begin match arg_dels with
           (* if no args, no delete instruction *)
           | [] ->
              trm_seq ~loc:t.loc
                (bodyl ++
                 [
                   trm_labelled return_label
                     (change_trm (trm_apps (trm_var f) arg_vals)
                        (trm_lit Lit_unit) t)
                 ]
                )
           | _ ->
              trm_seq ~annot:(Some Delete_instructions)
                ((trm_seq ~loc:t.loc
                    (bodyl ++
                     [
                       trm_labelled return_label
                         (change_trm (trm_apps (trm_var f) arg_vals)
                            (trm_lit Lit_unit) t)
                     ]
                    )
                 ) ::
                 arg_dels
                )
           end
        | _ ->
           trm_seq ~annot:(Some Delete_instructions) ~loc:t.loc
             ([
                trm_seq ~loc:t.loc
                  (arg_decls ++ (result_decl :: bodyl) ++
                   [
                     trm_labelled return_label
                       (change_trm
                          (trm_apps (trm_var f) arg_vals)
                          (trm_apps ~annot:(Some Heap_allocated)
                             (trm_unop Unop_get) [trm_var result])
                          t
                       )
                   ]
                  );
                trm_apps ~annot:(Some Heap_allocated) ~loc:t.loc ~is_instr:true
                  ~typ:(Some (typ_unit ())) (trm_unop (Unop_delete false))
                  [trm_var result]
               ] ++
               arg_dels
             )
      in
      (* clean up *)
      let t = group_decl_init t in
      let t = eliminate_goto_next t in
      let n = nb_goto return_label t in
      if n = 0 then delete_label return_label t else t
  in
  List.fold_left
    (fun t pl ->
      let p = List.flatten pl in
      let b = !Flags.verbose in
      Flags.verbose := false;
      let epl = resolve_path p t in
      Flags.verbose := b;
      match epl with
      | [] ->
         print_info t.loc "inline_fun_decl: no matching subterm for path %s\n"
           (string_of_path p);
         t
      | _ ->
         List.fold_left (apply_local_transformation apply_change) t epl
    )
    t
    inline_at

(*
  inline the definition pointed at by pl
  paths point at subterms in which all occurences will be replaced
  the empty path means all occurences will be replaced (default behaviour)
  optional argument to remove the declaration (not removed by default)
  assumption for function inlining: the function is used at most once per
  instruction
 *)
let inline_decl (clog : out_channel) ?(delete_decl : bool = false)
  ?(inline_at : path list list = [[]]) ?(fun_result : var = "res")
  ?(fun_return_label : label = "exit") (pl : path list) (t : trm) : trm =
  let p = List.flatten pl in
  let b = !Flags.verbose in
  Flags.verbose := false;
  let epl = resolve_path p t in
  Flags.verbose := b;
  match epl with
  | [dl] ->
     let t =
       let (t_def, _) = resolve_explicit_path dl t in
       let log : string =
         let loc : string =
           match t_def.loc with
           | None -> ""
           | Some (_, line) -> Printf.sprintf "at line %d " line
         in
         Printf.sprintf
           ("  - expression\n%s\n" ^^
            "    %sis a declaration\n"
           )
           (ast_to_string t_def) loc
       in
       write_log clog log;
       match t_def.desc with
       (* const variables *)
       | Trm_decl (Def_var ((x, _), dx)) ->
          let t_x = trm_var x in
          change_trm ~change_at:inline_at t_x dx t
       (*
         heap allocated variables
         note: an initialisation must be given
        *)
       | Trm_seq [{desc = Trm_decl (Def_var ((x, _), _)); _};
                  {desc = Trm_apps (_, [_; dx]); _}]
            when t_def.annot = Some Heap_allocated ->
          let t_x =
            trm_apps ~annot:(Some Heap_allocated) (trm_unop Unop_get)
              [trm_var x]
          in
          (*
             make sure x is not replaced in delete instructions by replacing it
             with a fresh variable
           *)
          let x' = fresh_in t x in
          let t =
            change_trm
              (trm_apps (trm_unop (Unop_delete false)) [trm_var x])
              (* do not forget annotations *)
              (trm_apps ~annot:(Some Heap_allocated) ~typ:(Some (typ_unit ()))
                 (trm_unop (Unop_delete false)) [trm_var x'])
              t
          in
          let t = change_trm ~change_at:inline_at t_x dx t in
          (* put back x instead of x' *)
          change_trm (trm_var x') (trm_var x) t
       (* typedef *)
       | Trm_decl (Def_typ (x, dx)) ->
          let ty_x = typ_var x in
          change_typ ~change_at:inline_at ty_x dx t
       (* fun decl *)
       | Trm_decl (Def_fun (f, tf, args, body)) ->
          let log : string =
            Printf.sprintf
              "  - function %s is used at most once per instruction\n"
              f
          in
          write_log clog log;
          inline_fun_decl ~inline_at fun_result fun_return_label f tf args body
            t
       | _ -> fail t.loc "inline_decl: expected a definition"
     in
     if delete_decl then remove_decl clog pl t else t
  | _ -> fail t.loc "inline_decl: the path must point at exactly 1 subterm"

let inline_seq (clog : out_channel) (pl : path list) (t : trm) : trm =
  let p = List.flatten pl in
  let b = !Flags.verbose in
  Flags.verbose := false;
  let epl = resolve_path p t in
  Flags.verbose := b;
  match epl with
  | [] ->
     print_info t.loc "inline_seq: no matching subterm";
     t
  | _ ->
     List.fold_left
       (fun t dl ->
         let log : string =
           let (t, _) = resolve_explicit_path dl t in
           let loc : string =
             match t.loc with
             | None -> ""
             | Some (_, line) -> Printf.sprintf "at line %d " line
           in
           Printf.sprintf
             ("  - expression\n%s\n" ^^
              "    %sis a seq inside another seq\n"
             )
             (ast_to_string t) loc
         in
         write_log clog log;
         let dl = List.rev dl in
         let n =
           match List.nth_opt dl 0 with
           | Some (Dir_nth n) -> n
           | _ -> fail t.loc "inline_seq: the path must point at a seq element"
         in
         apply_local_transformation
           (fun (t : trm) ->
             match t.desc with
             | Trm_seq tl ->
                let tl =
                  foldi
                    (fun i tl (t : trm) ->
                      if i <> n then tl ++ [t]
                      else
                        match t.desc with
                        | Trm_seq tl' -> tl ++ tl'
                        | _ -> fail t.loc "inline_seq: expected a seq"
                    )
                    []
                    tl
                in
                trm_seq ~annot:t.annot ~loc:t.loc ~add:t.add
                  ~attributes:t.attributes tl
             | _ -> fail t.loc "inline_seq: expected a seq container"
           )
           t
           (* remove the last direction to point at the seq *)
           (List.rev (List.tl dl))
       )
       t
       epl

let add_attribute (clog : out_channel) (a : attribute) (pl : path list)
  (t : trm) : trm =
  let p = List.flatten pl in
  let b = !Flags.verbose in
  Flags.verbose := false;
  let epl = resolve_path p t in
  Flags.verbose := b;
  match epl with
  | [] ->
     print_info t.loc "add_attribute: no matching subterm";
     t
  | _ ->
     List.fold_left
       (apply_local_transformation
          (fun t ->
            let log : string =
              let loc : string =
                match t.loc with
                | None -> ""
                | Some (_, line) -> Printf.sprintf "at line %d " line
              in
              Printf.sprintf
                ("  - expression\n%s\n" ^^
                 "    %sis a variable/type declaration\n"
                )
                (ast_to_string t) loc
            in
            write_log clog log;
            match t.desc with
            | Trm_decl (Def_var ((x, tx), init)) ->
               let ty_attributes = a :: tx.ty_attributes in
               {t with
                 desc = Trm_decl (Def_var ((x, {tx with ty_attributes}), init))}
            | Trm_decl (Def_typ (x, tx)) ->
               let ty_attributes = a :: tx.ty_attributes in
               {t with desc = Trm_decl (Def_typ (x, {tx with ty_attributes}))}
            | Trm_seq (t_decl :: tl) when t.annot = Some Heap_allocated ->
               begin match t_decl.desc with
               | Trm_decl (Def_var ((x, tx), init)) ->
                  begin match tx.ty_desc with
                  | Typ_ptr ty ->
                     let tx =
                       {tx with ty_desc =
                        Typ_ptr {ty with ty_attributes = a :: ty.ty_attributes}}
                     in
                     let t_decl =
                       {t_decl with desc = Trm_decl (Def_var ((x, tx), init))}
                     in
                     {t with desc = Trm_seq (t_decl :: tl)}
                  | _ -> assert false
                  end
               | _ -> assert false
               end
            | _ -> {t with attributes = a :: t.attributes}
          )
       )
       t epl
  
  (* 
    transforma loop of the shape
    optional label:
    for(int i = 0; i < N;i += D)
      body
    into a loop of the form 
    optional label:
    for (int c = 0; c < C; c++)
      for (int i = c*D; i < N; i += C*D)
      body
  *)

 
  let rec loop_coloring_aux (clog : out_channel) (c : var) (new_var : var)(t : trm) : trm =
    match t.desc with 
    (* The loop might be labelled, so keep the label *)
    | Trm_labelled (l, t_loop) ->
      trm_labelled l (loop_coloring_aux clog c new_var t_loop)
    
    | Trm_seq [t_loop; t_del] when t.annot = Some Delete_instructions ->
     let t_transformed = loop_coloring_aux clog c new_var t_loop in
     (* transformed loops are expected to declare their index *)
     begin match t_transformed.desc with
     | Trm_seq [{desc = Trm_for (init1, cond1, step1,
                                 {desc = Trm_seq [body1]; _}); _}; t_del1]
          when t_transformed.annot = Some Delete_instructions ->
        begin match body1.desc with
        | Trm_seq [{desc = Trm_for (init2, cond2, step2, body); _}; t_del2]
             when body1.annot = Some Delete_instructions ->
           (* if the index is used in body, then add delete instruction *)
           let i = deleted_var t_del in
           if not (is_used_var_in body i) then t_transformed
           else
             trm_seq ~annot:(Some Delete_instructions)
               [
                 trm_for init1 cond1 step1
                   (trm_seq
                      [trm_seq ~annot:(Some Delete_instructions)
                         [
                           trm_for init2 cond2 step2
                             (trm_seq ~annot:(Some Delete_instructions)
                                [body; t_del]);
                           t_del2
                         ]
                      ]
                   );
                 t_del1
               ]
        | _ -> fail body1.loc "loop_coloring_aux: expected inner loop"
        end
     | _ -> fail t_transformed.loc "loop_coloring_aux: expected outer loop"
     end 
  | Trm_for (init, cond, step, body) ->
      let log : string = 
        let loc : string = 
          match body.loc with
          | None -> ""
          | Some (_, line) -> Printf.sprintf "at line %d " line
        in 
        Printf.sprintf
         ("  - for (%s; %s; %s) is of the form\n" ^^
          "      for ([int] i = 0; i < N; i++)\n" ^^
          "  - expression\n%s\n" ^^
          "    %sis of the form\n" ^^
          "      {\n" ^^
          "        body\n" ^^
          "      }\n"
         )
         (ast_to_string init) (ast_to_string cond) (ast_to_string step)
         (ast_to_string body) loc    
      in
      write_log clog log;
      let index_i = for_loop_index t in
      let loop_size = for_loop_bound t in
      let block_size = trm_var c in
      let loop_step = for_loop_step t in 
      let log : string = 
        Printf.sprintf " - %s is divisible by %S\n" (ast_to_string loop_size) (ast_to_string block_size)
      in
      write_log clog log;
        
      let loop ?(top : bool = false) (index : var) (bound : trm) (body : trm) = 
        
        
        let start = match top with 
        | true -> trm_lit(Lit_int 0) 
        | false ->  
          match loop_step.desc with 
          | Trm_val(Val_lit(Lit_int 1)) -> trm_var new_var
          | _ -> trm_apps (trm_binop Binop_mul)
              [
                  trm_apps ~annot:(Some Heap_allocated) 
                      (trm_unop Unop_get) [trm_var new_var];
                    loop_step

              ]
              
          in 
          trm_seq ~annot:(Some Delete_instructions)
            [
              trm_for
                (*init *)
                (trm_seq ~annot:(Some Heap_allocated)
                  [
                    trm_decl (Def_var ((index, typ_ptr (typ_int())), trm_prim (Prim_new (typ_int ()))));
                    trm_set ~annot:(Some Initialisation_instruction)
                    (trm_var index) start
                  ]
                )
                (* cond *)
                (trm_apps (trm_binop Binop_lt)
                  [
                    trm_apps ~annot:(Some Heap_allocated)
                      (trm_unop Unop_get) [trm_var index];
                      bound
                  ]
                )
                (* step *)
                
                (if top then trm_apps (trm_unop Unop_inc) [trm_var index]
                else  match loop_step.desc with 
                  | Trm_val(Val_lit(Lit_int 1)) -> trm_set (trm_var index) ~annot:(Some App_and_set)
                    (trm_apps (trm_binop Binop_add)
                      [
                        trm_var index ;
                          
                        trm_var c
                      ])
                  | _ -> 
                    trm_set (trm_var index) ~annot:(Some App_and_set) (trm_apps (trm_binop Binop_add)
                      [
                        trm_var index;
                        trm_apps (trm_binop Binop_mul)
                             [
                               trm_apps ~annot:(Some Heap_allocated)
                                 (trm_unop Unop_get) [trm_var c];
                               loop_step
                             ]
                      
                      ])
                )
                  
                  
                

                (* body *)
                body;
                trm_apps ~annot:(Some Heap_allocated) ~typ:(Some (typ_unit ()))
                  (trm_unop (Unop_delete false)) [trm_var index]

            ]
          in loop ~top:true new_var (trm_var c) (trm_seq [loop ~top:false index_i loop_size body ]) 
      

  | _ -> fail t.loc "loop_coloring_aux: not a for loop, check the path "


let loop_coloring_aux (clog : out_channel)(c : var)(new_var : var) (t : trm) : trm = 
  let log : string = 
    let loc : string = 
      match t.loc with 
      | None -> ""
      | Some(_, line) -> Printf.sprintf "at line %d" line
    
    in
    Printf.sprintf
      ("  - expression\n%s\n" ^^

       "    %sis a (labelled) loop\n"
      )
      (ast_to_string t) loc 
    in
    write_log clog log;
    loop_coloring_aux clog c new_var t


let loop_coloring (clog : out_channel) (pl : path list) (c : var)(new_var : var)(t : trm) : trm = 
  let p = List.flatten pl in
  let b = !Flags.verbose in
  Flags.verbose := false;
  let epl = resolve_path p t in
  Flags.verbose := b;
  match epl with
  | [] ->
     print_info t.loc "loop_coloring: no matching subterm\n";
     t
  | _ ->
     List.fold_left
       (fun t dl ->
         apply_local_transformation (loop_coloring_aux clog c new_var) t dl)
       t
       epl

let rec loop_tile_aux (clog : out_channel)(b : var)(new_var : var) (t : trm) : trm =
  match t.desc with
  (* the loop might be labelled: keep the label *)
  | Trm_labelled (l, t_loop) ->
     trm_labelled l (loop_tile_aux clog b new_var t_loop)
  (*
    if the loop declares its own index, a seq with a delete instruction occurs
    in this case, put the delete instructions at the end of the inner loop if
    the index is still used
   *)
  | Trm_seq [t_loop; t_del] when t.annot = Some Delete_instructions ->
     let t_tiled = loop_tile_aux clog b new_var t_loop in
     (* tiled loops are expected to declare their index *)
     begin match t_tiled.desc with
     | Trm_seq [{desc = Trm_for (init1, cond1, step1,
                                 {desc = Trm_seq [body1]; _}); _}; t_del1]
          when t_tiled.annot = Some Delete_instructions ->
        begin match body1.desc with
        | Trm_seq [{desc = Trm_for (init2, cond2, step2, body); _}; t_del2]
             when body1.annot = Some Delete_instructions ->
           (* if the index is used in body, then add delete instruction *)
           let i = deleted_var t_del in
           if not (is_used_var_in body i) then t_tiled
           else
             trm_seq ~annot:(Some Delete_instructions)
               [
                 trm_for init1 cond1 step1
                   (trm_seq
                      [trm_seq ~annot:(Some Delete_instructions)
                         [
                           trm_for init2 cond2 step2
                             (trm_seq ~annot:(Some Delete_instructions)
                                [body; t_del]);
                           t_del2
                         ]
                      ]
                   );
                 t_del1
               ]
        | _ -> fail body1.loc "loop_tile_aux: expected inner loop"
        end
     | _ -> fail t_tiled.loc "loop_tile_aux: expected outer loop"
     end
  (* otherwise, just tile *)
  | Trm_for (init, cond, step, body) ->
     let log : string =
       let loc : string =
         match body.loc with
         | None -> ""
         | Some (_, line) -> Printf.sprintf "at line %d " line
       in
        Printf.sprintf
         ("  - for (%s; %s; %s) is of the form\n" ^^
          "      for ([int] x = 0; x < X; x++)\n" ^^
          "  - expression\n%s\n" ^^
          "    %sis of the form\n" ^^
          "      {\n" ^^
          "        body\n" ^^
          "      }\n"
         )
         (ast_to_string init) (ast_to_string cond) (ast_to_string step)
         (ast_to_string body) loc
    in
    write_log clog log;
    let index_x = for_loop_index t in 
    let loop_size = for_loop_bound t in 
    let block_size =  trm_var b in 
    let spec_bound = trm_apps (trm_var "min") 
          [
            loop_size;
            trm_apps (trm_binop Binop_add)
            [ 
             
              trm_var ("b" ^ index_x);
              trm_apps ~annot:(Some Heap_allocated)
                      (trm_unop Unop_get) [trm_var b]
            ]
          ]
    in
    let log : string = 
      Printf.sprintf "   -%s is divisible by %S\n" (ast_to_string loop_size) (ast_to_string block_size)
    in 
    write_log clog log;
  
    let loop ?(top : bool = false) (index : var) (bound : trm) (body : trm) =
        let start = match top with 
        | true -> trm_lit(Lit_int 0)
        | false -> trm_var( new_var)
        in 
        trm_seq ~annot:(Some Delete_instructions)
            [
              trm_for
                (* init *)
                (trm_seq ~annot:(Some Heap_allocated)
                   [
                     trm_decl (Def_var ((index, typ_ptr (typ_int ())),
                                        trm_prim (Prim_new (typ_int ()))));
                     trm_set ~annot:(Some Initialisation_instruction)
                       (trm_var index) start
                   ]
                )
                (* cond *)
                 (trm_apps (trm_binop Binop_lt)
                   [
                     trm_apps ~annot:(Some Heap_allocated)
                       (trm_unop Unop_get) [trm_var index];
                     bound
                   ]
                    )
                (* step *)
                (if not top then trm_apps (trm_unop Unop_inc) [trm_var index]
                else trm_set (trm_var index ) ~annot:(Some App_and_set)(trm_apps (trm_binop Binop_add)
                    [
                      trm_var index;
                      trm_apps ~annot:(Some Heap_allocated)
                      (trm_unop Unop_get) [trm_var b]

                    ]
                )
                
                )
                (* body *)
                body;
              trm_apps ~annot:(Some Heap_allocated) ~typ:(Some (typ_unit ()))
                (trm_unop (Unop_delete false)) [trm_var index]
            ]
        in
        loop ~top:true new_var loop_size (trm_seq [loop ~top:false index_x spec_bound body])
     | _ -> fail t.loc "loop_tile_aux: bad loop body"
    


let loop_tile_aux (clog : out_channel)(b : var)(new_var : var) (t : trm) : trm =
  let log : string =
    let loc : string =
      match t.loc with
      | None -> ""
      | Some (_, line) -> Printf.sprintf "at line %d " line
    in
    Printf.sprintf
      ("  - expression\n%s\n" ^^
       "    %sis a (labelled) loop\n"
      )
      (ast_to_string t) loc
  in
  write_log clog log;
  loop_tile_aux clog b new_var t

let loop_tile (clog : out_channel) (pl : path list)(tile_width : var)(new_var : var)(t : trm) : trm =
  let p = List.flatten pl in
  let b = !Flags.verbose in
  Flags.verbose := false;
  let epl = resolve_path p t in
  Flags.verbose := b;
  match epl with
  | [] ->
     print_info t.loc "loop_tile: no matching subterm\n";
     t
  | _ ->
     List.fold_left
       (fun t dl ->
         apply_local_transformation (loop_tile_aux clog tile_width new_var) t dl)
       t
       epl

let rec loop_swap_aux (clog : out_channel) (t : trm) : trm = 
  match t.desc with 
  (* the loop might be labelled: kepp the label *)
  | Trm_labelled (l, t_loop) -> 
    trm_labelled l (loop_swap_aux clog t_loop)
  | Trm_seq [t_loop; t_del] when t.annot = Some Delete_instructions ->
     let t_swaped = loop_swap_aux clog t_loop in
     (* swaped loops are expected to declare their index *)
     begin match t_swaped.desc with
     | Trm_seq [{desc = Trm_for (init1, cond1, step1,
                                 {desc = Trm_seq [body1]; _}); _}; t_del1]
          when t_swaped.annot = Some Delete_instructions ->
        begin match body1.desc with
        
        | Trm_seq [{desc = Trm_for (init2, cond2, step2, body); _}; t_del2]
             when body1.annot = Some Delete_instructions ->
           (* if the index is used in body, then add delete instruction *)
           let i = deleted_var t_del in
           if not (is_used_var_in body i) then t_swaped
           else
             trm_seq ~annot:(Some Delete_instructions)
               [
                 trm_for init1 cond1 step1
                   (trm_seq
                      [trm_seq ~annot:(Some Delete_instructions)
                         [
                           trm_for init2 cond2 step2
                             (trm_seq ~annot:(Some Delete_instructions)
                                [body; t_del]);
                           t_del2
                         ]
                      ]
                   );
                 t_del1
               ]
        | _ -> fail body1.loc "loop_swap_aux: expected inner loop"
        end
     | _ -> fail t_swaped.loc "loop_swap_aux: expected outer loop"
     end
  (* otherwise, just swap  *)
  | Trm_for (init1, cond1, step1,body1) ->
    let log : string = 
      let loc : string = 
        match body1.loc with 
        | None -> ""
        | Some (_, line) -> Printf.sprintf "at line %d " line
      in 
      Printf.sprintf
         ("  - for (%s; %s; %s) is of the form\n" ^^
          "      for ([int] x = 0; x < X; x++)\n" ^^
          "  - expression\n%s\n" ^^
          "    %sis of the form\n" ^^
          "      {\n" ^^
          "        body\n" ^^
          "      }\n"
         )
         (ast_to_string init1) (ast_to_string cond1) (ast_to_string step1)
         (ast_to_string body1) loc
      in
      write_log clog log;

      begin match body1.desc with 
     
      | Trm_seq ({desc = Trm_seq(f_loop :: _);_} :: _) ->
        begin match f_loop.desc with 
        | Trm_for(init2,cond2,step2,body2) -> 
          let log : string = 
            let loc : string = 
            match body1.loc with 
            | None -> ""
            | Some (_, line) -> Printf.sprintf "at line %d " line
          in 
          Printf.sprintf
          ("Inner looop " ^^
           "  - for (%s; %s; %s) is of the form\n" ^^
           "      for ([int] x = 0; x < X; x++)\n" ^^
           "  - expression\n%s\n" ^^
           "    %sis of the form\n" ^^
           "      {\n" ^^
           "        body\n" ^^
           "      }\n"
          )
          (ast_to_string init2) (ast_to_string cond2) (ast_to_string step2)
          (ast_to_string body2) loc
          in
          write_log clog log;
          let index1 = for_loop_index t in 
          let loop_size1 = for_loop_bound t in 
          let index_init1 = for_loop_init t in
          let index2 = for_loop_index f_loop in 
          let loop_size2 = for_loop_bound f_loop in
          let index_init2 = for_loop_init f_loop in
           
          let loop (index : var) (init : trm) (step : trm) (bound : trm) (body : trm) =
          trm_seq ~annot:(Some Delete_instructions)
            [
              trm_for
                (* init *)
                (trm_seq ~annot:(Some Heap_allocated)
                   [
                     trm_decl (Def_var ((index, typ_ptr (typ_int ())),
                                        trm_prim (Prim_new (typ_int ()))));
                     trm_set ~annot:(Some Initialisation_instruction)
                       (trm_var index) (init)
                   ]
                )
                (* cond *)
                (trm_apps (trm_binop Binop_lt)
                   [
                     trm_apps ~annot:(Some Heap_allocated)
                       (trm_unop Unop_get) [trm_var index];
                     bound
                   ]
                )
                (* step *)
                (step)
                (* body *)
                body;
              trm_apps ~annot:(Some Heap_allocated) ~typ:(Some (typ_unit ()))
                (trm_unop (Unop_delete false)) [trm_var index]
            ]
        in
        loop index2 index_init2 step2 loop_size2 (trm_seq [loop index1 index_init1 step1 loop_size1 body2])
        | _ -> fail t.loc "loop_swap_aux: inner_loop was not matched"
        end
      | _ -> fail t.loc "loop_swap_aux; expected inner loop"
      end

      
  | _ -> fail t.loc "loop_swap_aux; bad loop body"

  

let loop_swap_aux (clog : out_channel) (t : trm) : trm = 
  let log : string =
    let loc : string =
      match t.loc with 
      | None -> ""
      | Some (_, line) -> Printf.sprintf "at line %d " line
    in 
    Printf.sprintf
      (" - expression \n%s\n" ^^
       "  %sis a (labelled) loop\n"
      )
      (ast_to_string t) loc
  in
  write_log clog log;
  loop_swap_aux clog t


let loop_swap (clog : out_channel) (pl : path list)(t : trm) : trm =
  let p = List.flatten pl in
  let b = !Flags.verbose in
  Flags.verbose := false;
  let epl = resolve_path p t in
  Flags.verbose := b;
  match epl with
  | [] ->
     print_info t.loc "loop_swap: no matching subterm\n";
     t
  | _ ->
     List.fold_left
       (fun t dl ->
         apply_local_transformation (loop_swap_aux clog) t dl)
       t
       epl
(* get_loop_nest_indices -- currently omiting the last one 
   
*)
(* for a { for b {} {  for j {}   ;  for k {} } } -- >  a::b::[]
  the function should check that it is a loop nest :
      aux t = 
         if t is a for (i , body) then  i::(aux body)
         if t is a seq[(for(..)) as t1] and nothing else then  aux t1
         else []
         
*)
let rec get_loop_nest_indices (t : trm) : 'a list = 
    match t.desc with 
    | Trm_labelled (_, t_loop) -> get_loop_nest_indices t_loop
    | Trm_seq [t_loop;_] -> get_loop_nest_indices t_loop
    | Trm_for (_,_,_,body) -> 
      let loop_index = for_loop_index t in 
      begin match body.desc with 
      | Trm_seq ({desc = Trm_seq (f_loop :: _);_} :: _) -> 
        loop_index :: get_loop_nest_indices f_loop
      | _ -> loop_index :: []
      end
    | _ -> []


let move_loop_before_aux (clog : out_channel) (loop_index : var) (t : trm) : trm =
    let log : string = 
      let loc : string = 
        match t.loc with 
        | None -> ""
        | Some (_, line) -> Printf.sprintf "at line %d " line
      in 
      Printf.sprintf
      ("  - expression\n%s\n" ^^
          "    %sis a (labelled) loop\n"
      )
      (ast_to_string t) loc
      in
      write_log clog log;
      (* Get the path from the outer loop to the one we want to swap with 

              
      let path_list = List.rev (get_path clog t) in 

        *)
      let path_list = List.rev (get_loop_nest_indices t)  in      
      (* do a list rev at the end of get_loop_nest_indices 
         let rec chop_list_before x xs =
            | [] -> error "did not find x"
            | y::tl -> if y = x then [] else y:: chop_list_before x tl      *)
      let rec clean_path (xl : 'a list) : 'a list = match xl with 
        | [] -> []
        | hd :: tl -> 
          if hd = loop_index then tl
          else clean_path tl
      in 
      let _check_last = List.mem loop_index path_list in 
      (*
      let path_list = if not check_last then path_list
          else clean_path path_list 
      *)
      let path_list = clean_path path_list 
      in
      (* List.fold_right (fun i acc  -> loop swap t i) path_list acc t 
         --checkout the documentation of fold_right *)
      let rec multi_swap (xl : 'a list) (t : trm) : trm = match xl with 
      | [] -> t
      | hd :: tl -> 
        let pl = [cFor ~init:[cVarDef ~name:hd ()] ()] in
        let t = loop_swap clog pl t in 
        multi_swap tl t
     in
     multi_swap path_list t 

let move_loop_before (clog : out_channel) (pl : path list)(loop_index : var) (t : trm) : trm =
  let p = List.flatten pl in 
  let b = !Flags.verbose in
  Flags.verbose := false;
  let epl = resolve_path p t in 
  Flags.verbose := b;
  match epl with 
  | [] -> 
    print_info t.loc "move_loop_before: no matching subterm";
    t
  | _ ->
    List.fold_left
      (fun t dl ->
        apply_local_transformation (move_loop_before_aux clog loop_index) t dl)
      t
      epl

let move_loop_after_aux (clog : out_channel) (loop_index : var) (t : trm) : trm = 
  let log : string =
    let loc : string =
      match t.loc with 
      | None -> ""
      | Some (_, line) -> Printf.sprintf "at line %d " line
    in
    Printf.sprintf
    ("  - expression\n%s\n" ^^
          "    %sis a (labelled) loop\n"
    )
    (ast_to_string t) loc
    in
    write_log clog log;
    let path_list = get_loop_nest_indices t in
    let rec clean_path (xl : 'a list) : 'a list = match xl with 
      | [] -> []
      | hd :: tl -> 
        if hd = loop_index then tl
        else clean_path tl  
      in
    let l_index = List.hd path_list in
    let _check_last = List.mem loop_index path_list in 
    let path_list = if false then path_list
      else clean_path (List.rev path_list) 
      in
    let path_length = List.length path_list in
    (*if (path_list = []) then error ---try to check the error in case  move_before "i" "i" *)
    (* List.fold_left (fun _i acc -> swap l_index) (List.tl path_list) *)
    let rec multi_swap (count : int) (t : trm) : trm = match count with 
      | 0 ->  t  
      | _ -> let pl = [cFor ~init:[cVarDef ~name:l_index ()] ()] in 
           let t = loop_swap clog pl t in 
           multi_swap (count-1) t
      in
    multi_swap path_length t 

let move_loop_after (clog : out_channel) (pl : path list)(loop_index : var) (t : trm) : trm =
  let p = List.flatten pl in 
  let b = !Flags.verbose in
  Flags.verbose := false;
  let epl = resolve_path p t in 
  Flags.verbose := b;
  match epl with 
  | [] -> 
    print_info t.loc "move_loop_before: no matching subterm";
    t
  | _ ->
    List.fold_left
      (fun t dl ->
        apply_local_transformation (move_loop_after_aux clog loop_index) t dl)
      t
      epl

let move_loop (clog : out_channel)  ?(move_before : string = "") ?(move_after : string = "") (loop_index : string) (t : trm) : trm = 
  let log : string = 
      let loc : string = 
        match t.loc with 
        | None -> ""
        | Some (_, line) -> Printf.sprintf "at line %d " line
      in Printf.sprintf
          ("  - expression\n%s\n" ^^
          "    %sis a struct type\n"
          )
      (ast_to_string t) loc 
    in
  write_log clog log;
  match move_before, move_after with 
  | "",_ -> move_loop_after clog [cFor ~name:loop_index ()] move_after t 
  | _,"" -> move_loop_before clog [cFor ~name:move_before ()] loop_index t
  | _,_ -> fail t.loc "move_loop: only one of move_before or move_after should be specified"
  


let find_keys value m = 
  Field_map.fold(fun k v acc -> if v = value then k :: acc else acc) m []

(* A function rename all th elements of a list *)
let rec apply_labels vl pl = match pl with 
| [] -> []
| hd :: tl -> let y = List.map (fun x -> hd ^ "_" ^x) vl in y :: apply_labels vl tl


let add_key key value m = Field_map.add key value m 


let rec add_keys (lk : var list) (lv : typ list) m  = match (lk ,lv) with 
| [],[] -> m
| _ :: _, [] -> m
| [] , _ :: _ -> m
| hd :: tl, hd1 :: tl1 -> let m = add_key hd hd1 m in add_keys tl tl1 m;;

let rec add_keys_to_map lv llk m = match llk with 
| [] -> m 
| hd :: tl -> let m = add_keys  hd lv m in add_keys_to_map lv tl m

(*
let record_get_typed_fields (fields_list, fields_map) =
                list (string * typ) : list =
                 List.combine fields_list (get_values fields_list fields_map
*)


let rec insert_before x local_l list = match list with 
| [] -> []
| hd :: tl -> if hd = x then local_l @ hd :: tl else hd :: (insert_before x local_l tl)

let rec insert_list keys_list temp_field_list field_list1 = match keys_list with 
| [] -> field_list1
| hd :: tl -> let field_list1 = insert_before hd (List.hd temp_field_list) field_list1 in insert_list tl (List.tl temp_field_list ) field_list1


let inline_struct_aux (clog : out_channel) ?(struct_fields : fields = []) (t1 : trm) (t : trm) : trm =
  let log : string = 
    let loc : string = 
      match t.loc with
      | None -> ""
      | Some (_, line) -> Printf.sprintf "at line %d " line
    in Printf.sprintf
      ("  - expression\n%s\n" ^^
      "    %sis a struct type\n"
      )
    (ast_to_string t) loc 
    in 
    write_log clog log;
    begin match t1.desc with 
      | Trm_decl (Def_typ(_,dx)) -> 
        let field_list, field_map =
          match dx.ty_desc with 
            | Typ_struct (l,m,_) -> l,m
            | _ -> fail t.loc "inline_struct_aux: The type shoudl be a typedef struct"
        in  
        begin match t.desc with 
        | Trm_decl (Def_typ(x1,dx1)) -> 
            let field_list1, field_map1,name = 
              match dx1.ty_desc with
              | Typ_struct(l,m,n) -> l,m,n
              |_ -> fail t.loc "inline_struct_aux: the type should be a typedef struct"
            in 
          
          (* If the list of fields is given then do nothing otherwise find all occurrences of typedef first struct*)
          let keys_list = if struct_fields = [] then Field_map.fold(fun k v acc -> if v = (typ_var x1) then k :: acc else acc) field_map1 []

            else struct_fields 
            in
          
          (* Apply the labels *)
          let temp_field_list = apply_labels field_list keys_list in 
          
          (* The key values from the first struct *)
          let values = List.map (fun x -> Field_map.find x field_map) field_list in 

          (* Add the new keys with their values to the second  struct field_map *)
          let field_map1 = add_keys_to_map values temp_field_list field_map1 in 
          

          let field_list1 = insert_list keys_list temp_field_list field_list1 in 
          
          
          let _field_map1 = List.fold_left (fun mapPrev key -> Field_map.remove key mapPrev) field_map1 keys_list in
          
          let field_list1 = list_remove_set  keys_list field_list1 in 

          (* do removal at the end_*)
          (* (m',l') = remove_keys_from_list_and_map xs  (m,l) *)
          (*   insert_bindings_in_list_and_map  field_before_x  new_bindings (m,l) *)
          (*  List.fold_left (fun (x,t) acc -> Fmap.add x t acc) new_bindings m *)

          (* record_get_typed_fields (fields_list, fields_map) =
                list (string * typ) : list =
                 List.combine fields_list (get_values fields_list fields_map) *)

        (* Solution better: 
              let fmap_of_list xts =  List.fold_left (fun (x,t) acc -> Fmap.add x t acc) new_bindings Fmap.empty 
              typ_struct typ_fields_list = Typ_struct (List.keys typ_fields_list, Fmap.of_list typ_fields_list)
        *)
          (* List documentation 
          let prefix = here it is "pos"
          let xts = (record_get_typed_fields .. ) in
          let yts = List.map (fun (x,t) -> (prefix ^ "_" ^ x, t)) xts in
          insert_bindings_in_list prefix yts target_bindings

          let insert_bindings_in_list prefix yts xts =
            match xts with
            | [] -> error 
            | x::xts' -> if x = prefix then x::(yts@xts') else x::(insert_bindings_in_list prefix yts xts')
            

           List;fold left right f
           List.fold_left2   (List.combine / split)  *)

          trm_decl (Def_typ (x1,typ_struct field_list1 field_map1 name))
        
        | _ -> fail t.loc "inline_struct_aux: expected a definiton"
        end
      | _ -> fail t.loc " inline_struct_aux: expected a definiton"
      end




let change_struct_access  (x : typvar) (t : trm) : trm = 
  let rec aux (global_trm : trm) (t : trm) : trm = 
    match t.desc with 
    | Trm_apps (f, [base]) ->
       begin match f.desc with 
        | Trm_val (Val_prim (Prim_unop (Unop_struct_access y)))
          | Trm_val (Val_prim (Prim_unop (Unop_struct_get y))) ->
        begin match base.desc with 
        | Trm_apps (f',[base']) ->
          begin match f'.desc with 
          | Trm_val (Val_prim (Prim_unop (Unop_struct_access pos)))
            | Trm_val (Val_prim (Prim_unop (Unop_struct_get pos))) when pos = x -> 
            let new_var = pos ^"_"^ y
             in
            let new_f = {f' with desc = Trm_val(Val_prim (Prim_unop (Unop_struct_access new_var)))}
            in
            trm_apps ~annot:t.annot ~loc:t.loc ~is_instr:t.is_instr
                     ~add:t.add ~typ:t.typ new_f [base']
            | _ -> trm_map (aux global_trm) t
            end
          | _ -> trm_map (aux global_trm) t
          end
        | _ -> trm_map (aux global_trm) t
        end
    
      (* other cases: recursive call *)
    | _ -> trm_map (aux global_trm) t
in 
aux t t

let rec inline_sublist_in_list sublist i xs = match xs with 
| [] -> []
| h :: t -> if i = 0 then sublist @ t else h :: inline_sublist_in_list sublist (i-1) t

let rec copy_list l = match l with 
| [] -> []
| hd :: tl -> hd :: hd :: copy_list tl;;


(* Get the index for a given field of struct inside its list of fields *)
let get_pos (x : typvar) (t : trm) : int = 
  begin match t.desc with 
  | Trm_decl (Def_typ(_,dx)) -> 
       let field_list1 = 
          match dx.ty_desc with
          | Typ_struct(l,_,_) -> l
          |_ -> fail t.loc "get_pos: the type should be a typedef struct"
        in

        let rec find x lst = 
        match lst with 
        | [] -> raise (Failure "Not Found")
        | hd :: tl -> if hd = x then 0 else 1 + find x tl
        in 
        find x field_list1 
    | _ -> fail t.loc "get_pos_and_element: expected a struct type"
    end




let change_struct_initialization (_clog : out_channel) (struct_name : typvar) (base_struct_name : typvar) (x : typvar) (t :trm) : trm = 
  let base_struct_path = [cType ~name:base_struct_name()] in 
  let epl_of_base_struct = resolve_path (List.flatten base_struct_path) t in 
  let base_struct_term = match epl_of_base_struct with 
    | [dl] -> let (t_def,_) = resolve_explicit_path dl t in t_def 
    | _ -> fail t.loc "change_struct_initialization: expected a typedef struct"
  in 
  let struct_path = [cType ~name:struct_name ()] in 
  let epl_of_struct = resolve_path (List.flatten struct_path) t in 
  let struct_term = match epl_of_struct with 
  | [dl] -> 
    let (t_def,_) = resolve_explicit_path dl t in t_def 
  | _ -> fail t.loc "change_struct_initialization: expected a typedef struct"
  in 
  
  let pos = get_pos x struct_term in 
  let rec aux (global_trm : trm) (t : trm) = 
    match t.desc with 
    | Trm_struct term_list -> 

      begin match t.typ with 
      | Some{ ty_desc = Typ_var y;_} when y = struct_name -> 
        
        let el = List.nth term_list pos in 
        (*let show_trm t = Print_ast.print_ast ~only_desc:true stdout t;
        Dev.ml    Dev.trm t TODO 
        *)
        Print_ast.print_ast ~only_desc:true stdout el;
        begin match el.desc with 
        | Trm_struct inner_term_list -> trm_struct (inline_sublist_in_list inner_term_list pos term_list)
          
        | Trm_apps(_,[body]) -> 
          
          begin match body.desc with 
          | Trm_var _p ->  (*trm_struct (List.rev term_list)*) 
              let field_list = 
              match base_struct_term.desc with 
                | Trm_decl(Def_typ(_,dx)) ->
                  begin match dx.ty_desc with 
                    | Typ_struct (fl,_,_) -> fl 
                    | _ -> fail t.loc "change_struct_initializaition: expected a struct"
                  end
                | _ -> fail t.loc "change_struct_initialization: expected a definition"
              in 
              let field_list = List.map (fun el -> trm_var (_p ^ "." ^ el)) field_list
              in trm_struct (inline_sublist_in_list field_list pos term_list)
          | _ -> fail t.loc "change_struct_initialization: expected either a record or a variables"
          end
        | _ -> fail t.loc "change_struct_initialization: expected either a record or a variables"
        end
      | _ -> trm_map (aux global_trm) t
      end
    | _ -> trm_map (aux global_trm) t
  in 
  aux t t

        

let inline_struct (clog : out_channel)  ?(struct_fields : fields = []) (name : string)(t : trm) : trm =
  let field_name = List.hd struct_fields in 

  let struct_term_path  = [cType ~name:name ()] in 
  let p_of_struct_term = List.flatten struct_term_path in
  let epl_of_struct_term = resolve_path p_of_struct_term t in 
  let struct_term = match epl_of_struct_term with 
    | [dl] -> 
      let(t_def,_) = resolve_explicit_path dl t in t_def
    | _ -> fail t.loc "inline_struct: expected a typedef struct"
    in 
  (* Get the type of the field_name by going through the field_map of struct obj *)
  let inner_struct_name = 
  match struct_term.desc with 
  | Trm_decl (Def_typ(_,dx)) ->
    let field_map = 
      match dx.ty_desc with
      | Typ_struct (_,m,_) -> m
      | _ -> fail t.loc "inline_struct: the type should be a typedef struct"
    in 
    let field_map_typ = Field_map.find field_name field_map in 
    begin match field_map_typ.ty_desc with 
    | Typ_var y -> y
    | _ -> fail t.loc "inline_struct: expeted a typ var as the value of a key  "
    end
  | _ -> fail t.loc "inline_one_struct: expected a definition"
  in 
  
    
    let  pl_temp = [cType ~name:inner_struct_name()]  in
    let p_temp = List.flatten pl_temp in
    let epl_temp = resolve_path p_temp t in 
    let t1 = 
    match epl_temp with
    | [dl] -> 
      let (t_def,_) = resolve_explicit_path dl t in t_def
    | _ -> fail t.loc "inline_struct: expected a typedef struct"
    in
    
    
   let t =  List.fold_left (fun acc_t x -> change_struct_access x acc_t) t struct_fields
    in
    
    let t = List.fold_left (fun acc_t x -> change_struct_initialization  clog  name inner_struct_name x acc_t ) t struct_fields
    in
    
    let p = List.flatten struct_term_path in
    let b = !Flags.verbose in
    Flags.verbose := false;
    let epl = resolve_path p t in
    Flags.verbose := b;
    match epl with 
    | [] -> 
      print_info t.loc "inline_struct: no matching subterm";
      t
    | _ -> 
      List.fold_left 
        (fun t dl -> 
          apply_local_transformation (inline_struct_aux clog ~struct_fields t1) t dl)
        t
        epl 


  let _make_explicit_record_assignment (_clog : out_channel) (x : typvar) (t : trm) : trm =
    let struct_def_path = [cType ~name:x ()] in 
    let epl_of_struct_def_path = resolve_path (List.flatten struct_def_path) t in 
    let struct_def_term = match epl_of_struct_def_path with
    | [dl] -> let (t_def,_) = resolve_explicit_path dl t in t_def 
    | _ -> fail t.loc "make_explicit_record_assigment: expected a typedef struct"
    in 
    let field_list = 
    
    match struct_def_term.desc with
    | Trm_decl (Def_typ (_,dx)) -> 
      begin match dx.ty_desc with 
      | Typ_struct (fl,_,_) -> List.rev fl
      | _ -> fail t.loc "make_explicit_record_assigment: the type should be a struct" 
      end
    | _ -> fail t.loc "make_explicit_record_assigment: expected a definition"
    
    in
    let rec aux (global_trm : trm) (t : trm) : trm = 
      match t.desc with 
      | Trm_apps (f,[lt;rt]) ->
        begin match f.desc with 
        | Trm_var "overloaded=" ->
          begin match lt.desc with 
          | Trm_apps(f1,[lbase]) -> 
            begin match rt.desc with 
            Trm_apps (f2,[rbase]) -> 
                  let tl = List.map(fun sf -> 
                    let new_f = {f with desc = Trm_val(Val_prim (Prim_unop (Unop_struct_get sf)))}
                    in trm_apps ~annot:t.annot ~loc:t.loc ~is_instr:t.is_instr ~add:t.add ~typ:t.typ 
                    f [trm_apps ~annot:(Some Access) f1  [trm_apps new_f [lbase]];trm_apps ~annot:(Some Access) f2 [trm_apps new_f [rbase]]]
                  ) field_list in
                  trm_seq ~annot:(Some No_braces) tl
                  
            | _ -> fail t.loc "make_explicit_record_assigment: right term was not matched"
            end
          | _ -> fail t.loc "make_exmplicit_record_assigment: left term was not matched"
          end
        | _ -> trm_map (aux global_trm) t
        end
      | _ -> trm_map (aux global_trm) t 
    in aux t t


let make_explicit_record_assigment_aux (clog : out_channel) (field_list : fields) (t : trm) : trm = 
  let log : string = 
    let loc : string = 
     match t.loc with 
     | None -> ""
     | Some (_, line) -> Printf.sprintf "at line %d" line
    in Printf.sprintf
    (" -expression\n%s\n" ^^
    "   %sis an assigment\n"
    )
    (ast_to_string t) loc
    in write_log clog log;
    Print_ast.print_ast ~only_desc:true stdout t;
    begin match t.desc with
    | Trm_seq [t_decl;t_assign] -> 
      
      begin match t_assign.desc with 
      | Trm_apps(f,[lt;rt]) ->
        begin match f.desc with 
        | Trm_val( Val_prim (Prim_binop Binop_set )) -> 
          begin match rt.desc with 
          | Trm_apps (f',[rbase]) ->
            let tl = List.map(fun sf -> 
            let new_f = { f with desc = Trm_val(Val_prim (Prim_unop (Unop_struct_get sf)))}
            in trm_apps ~annot:t.annot ~loc:t.loc ~is_instr:t.is_instr ~add:t.add ~typ:t.typ 
            f [trm_apps ~annot:(Some Access) new_f [lt] ; trm_apps ~annot:(Some Access) f' [trm_apps new_f [rbase]]]

            ) field_list in 
            let tl = t_decl :: tl in 
            trm_seq ~annot:(Some No_braces) tl 
          | _ -> fail t.loc "make_explicit_record_assignment_aux: right term was not matched for the un-detached assigment case"
          end
        | _ -> fail t.loc "make_explicit_record_assignment_aux: function description was not matched "
        end
      | _ -> fail t.loc "make_explicit_record_assignment_aux: left term was not matched for the un-detached assignment case"
      end    
    
    | Trm_apps (f,[lt;rt]) ->

      begin match f.desc with 
      | Trm_var "overloaded=" ->
        begin match lt.desc with 
        | Trm_apps(f1,[lbase]) -> 
          begin match rt.desc with 
          Trm_apps (f2,[rbase]) -> 
            let tl = List.map(fun sf -> 
            let new_f = {f with desc = Trm_val(Val_prim (Prim_unop (Unop_struct_get sf)))}
            in trm_apps ~annot:t.annot ~loc:t.loc ~is_instr:t.is_instr ~add:t.add ~typ:t.typ 
            f [trm_apps ~annot:(Some Access) f1  [trm_apps new_f [lbase]];trm_apps ~annot:(Some Access) f2 [trm_apps new_f [rbase]]]
              ) field_list in 
              trm_seq ~annot:(Some No_braces) tl
                  
          | _ -> fail t.loc "make_explicit_record_assigment_aux: right term was not matched"
          end
        | _ -> fail t.loc "make_explicit_record_assigment_aux: left term was not matched"
        end
      | _ -> fail t.loc "make_explicit_record_assigment_aux: function name was not matched" 
      end

    | _ -> fail t.loc "make_explicit_record_assigment_aux: No variable declaration or assigment was matched"
    end




let make_explicit_record_assigment (clog : out_channel) ?(struct_name : string = "") (pl : path list) (t : trm) : trm = 
  let struct_def_path = [cType ~name:struct_name ()] in 
  let epl_of_struct_def_path = resolve_path (List.flatten struct_def_path) t in 
  let struct_def_term = match epl_of_struct_def_path with
  | [dl] -> let (t_def,_) = resolve_explicit_path dl t in t_def 
  | _ -> fail t.loc "make_explicit_record_assigment: expected a typedef struct"
  in 
  let field_list = 
    
  match struct_def_term.desc with
  | Trm_decl (Def_typ (_,dx)) -> 
    begin match dx.ty_desc with 
    | Typ_struct (fl,_,_) -> List.rev fl
    | _ -> fail t.loc "make_explicit_record_assigment: the type should be a struct" 
    end
  | _ -> fail t.loc "make_explicit_record_assigment: expected a definition"
  in   
  let p = List.flatten pl in 
  let b = !Flags.verbose in
  Flags.verbose := false; 
  let epl = resolve_path p t in 
  Flags.verbose := b;
  match epl with 
  | [] -> 
    print_info t.loc "make_explicit_record_assigment: no matching subterm";
    t
  |_ -> 
    List.fold_left 
      (fun t dl -> 
        apply_local_transformation (make_explicit_record_assigment_aux clog field_list) t dl)
        t 
        epl


