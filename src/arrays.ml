open Ast
open Target
open Ast_to_c
open Generic
open Arrays_core
open Output



let array_to_variables (tg : target) (new_vars : var list) : unit =
  apply_to_targets tg (fun p t ->
    (* TODO: Check together with Arthur, this particular solution *)
    let epl = resolve_target tg t in 
    let array_variable = 
    begin match epl with 
    | [dl] -> let (t_def, _) =  resolve_path dl t in 
      decl_name t_def 
    | _ -> fail t.loc "expected only one explicit path"
    end
    in
    let t = inline_array_access array_variable new_vars t in
    Arrays_core.array_to_variables new_vars p t)

(* TODO: Finish splititng all function into core and basics for this module *)

(*
  array tiling: transforms an array t[n] into a matrix t[n/b][b] for a fixed
  block size b
  name the block's type as block_name
  arguments:
    - x: type variable representing the array type
      x may be an alias for ty[n] or for ty* where ty is the type of the array
      elements
    - b: block size
  assumptions:
    - if x is ty*, each array of type x is allocated through a custom function:
      x a = my_alloc(nb_elements, size_element)
    - x is not used in function definitions, but only in var declarations
    - for now: in any case, the number of elements is divisible by b
 *)
let tile_array (clog : out_channel) (name : var -> var) (block_name : typvar)
  (b : trm) (x : typvar) (t : trm) : trm =
  (*
    changes:
      - replace the definition of x with:
        + (ty[b])* if x is ty*
        + ty[n/b][b] if x is ty[n]
        in both cases: define ty[b] as block_name
      - add a copy of each function taking an argument of type x and replace the
        function calls with these copies
      - replace array accesses a[i] for a of type x with a[i/b][i%b]
   *)
  let ilsm = functions_with_arg_type x t in
  (* first add copies of the functions *)
  let t = insert_fun_copies name ilsm x t in
  (* then replace function calls *)
  let t = replace_fun_names name ilsm x t in
  (* finally adapt the declaration and accesses *)
  let base_type =
    match aliased_type x t with
    | None -> fail t.loc "tile_array: unable to find array type"
    | Some ty ->
       let log : string =
         Printf.sprintf
           ("  - type\n%s\n" ^^
            "    represents a one-dimensional array type or a pointer type\n"
           )
           (typ_to_string ty)
       in
       write_log clog log;
       begin match ty.ty_desc with
       | Typ_ptr ty ->
          let log : string =
            Printf.sprintf
              ("  - each array of type %s is declared using the following " ^^
                 "pattern:\n" ^^
               "      %s array_name = my_alloc(nb_elements, size_element)\n" ^^
               "  where nb_elements is divisible by %s\n"
              )
              x x (ast_to_string b)
          in
          write_log clog log;
          ty
       | Typ_array (ty, s) ->
          let log : string =
            let n : string =
              match s with
              | Undefined ->  fail t.loc "tile_array: array size must be provided"
              | Const n -> string_of_int n
              | Trm t -> ast_to_string t
            in
            Printf.sprintf "  - size %s is divisible by %s\n" n
              (ast_to_string b)
          in
          write_log clog log;
          ty
       | _ -> fail t.loc "tile_array: expected array or pointer type"
       end
  in
  clean_up_no_brace_seq (tile_array_core base_type block_name b x t)


(*  t[i][k]  t:x   then swap dimentions for t

  transformation to swap the two first dimensions of an array
  name is used to name function copies
  assumption: x is a type variable that represents a multidimensional array type
  with >= 2 dimensions
  all variables of type x will be swapped
  assumption: x is not used in fun declarations
*)

(*
  swap the dimensions in the declaration of x
  swap all the accesses to arrays of type x
 *)

let array_swap (name : var -> var) (x : typvar) (t : trm) : trm =
  (*
    3 things to change in t:
    - dimensions in the declaration of x
    - add a copy of each function taking an argument of type x and replace the
    function calls with these copies
    - array accesses on variables of type x
      includes initialisation (through loop)
      --> no initialisation by list for multidimensional arrays
   *)
  let ilsm = functions_with_arg_type x t in
  (* first add copies of the functions *)
  let t = insert_fun_copies name ilsm x t in
  (* then replace function calls *)
  let t = replace_fun_names name ilsm x t in
  (* finally adapt the declaration and accesses *)
  Arrays_core.array_swap_aux x t


(*
  Array of Structures to Structure of Arrays:²²²
  if s is struct {t1 field1; …; tm fieldm} and x is s[n], transforms x into s'
  where s' is struct {t1 field1[n]; …; tm fieldm[n]}
  arguments:
    - x: type variable representing the array type
    - name: to name function copies
  assumptions:
    - x is not used in function definitions, but only in var declarations
 *)

let swap_accesses (clog : out_channel) (x : typvar) (t : trm) : trm =
  let rec aux (global_trm : trm) (t :trm) : trm =
    match t.desc with
    (* declarations *)
    | Trm_typedef d ->
       begin match d with
       (* we have to change the declaration of x *)
       | Typedef_abbrev (y, ty) when y = x ->
        let log : string =
          Printf.sprintf
           ("  - type\n%s\n" ^^
            "    is an array type\n"
           )
           (typ_to_string ty)
        in
        write_log clog log;
          (*
            ty must be an array type over a struct type denoted by a type var
           *)
          begin match ty.ty_desc with
          | Typ_array ({ty_desc = Typ_var (y, _); _}, s) ->
             begin match aliased_type y global_trm with
             | None ->
                fail t.loc "swap_accesses: cannot find underlying struct type"
             | Some ty' ->
                let log : string =
                  Printf.sprintf
                    ("  - type\n%s\n" ^^
                       "    is a struct type\n"
                    )
                    (typ_to_string ty')
                in
                write_log clog log;
                begin match ty'.ty_desc with
                | Typ_struct (l,m, n) ->
                   let m =
                     Field_map.map
                       (fun ty'' ->
                         typ_array ~ty_attributes:ty.ty_attributes ty'' s) m
                   in
                   trm_typedef (Typedef_abbrev(x, typ_struct l m n))
                | _ ->
                   fail t.loc "swap_accesses: expected underlying struct type"
                end
             end
          | _ -> fail t.loc "swap_accesses: expected array type declaration"
          end
       (*
         other cases: type declarations (not x), fun declarations, var
         declarations (not of type x)
         arrays of type x are heap allocated
        *)
       | _ -> trm_map (aux global_trm) t
       end
    (* accesses: y[i].f becomes (y.f)[i] *)
    | Trm_apps (f, [base]) ->
       begin match f.desc with
       | Trm_val (Val_prim (Prim_unop (Unop_struct_access _)))
         | Trm_val (Val_prim (Prim_unop (Unop_struct_get _))) ->
          begin match base.desc with
          | Trm_apps (f', [base'; index]) ->
             begin match f'.desc with
             | Trm_val (Val_prim (Prim_binop Binop_array_access))
               | Trm_val (Val_prim (Prim_binop Binop_array_get)) ->
                (*
                  swap accesses only if the type of base' is x (or x* in case of
                  an access on a heap allocated variable)
                 *)
                begin match base'.typ with
                | Some {ty_desc = Typ_var (y, _); _} when y = x ->
                   (* x might appear both in index and in base' *)
                   let base' = aux global_trm base' in
                   let index = aux global_trm index in
                   (* keep outer annotations *)
                   trm_apps ~annot:t.annot ~loc:t.loc ~is_statement:t.is_statement
                     ~add:t.add ~typ:t.typ f' [trm_apps f [base']; index]
                | Some {ty_desc = Typ_ptr {ty_desc = Typ_var (y, _); _}; _}
                     when y = x ->
                   (* x might appear both in index and in base' *)
                   let base' = aux global_trm base' in
                   let index = aux global_trm index in
                   (* keep outer annotations *)
                   trm_apps ~annot:t.annot ~loc:t.loc ~is_statement:t.is_statement
                     ~add:t.add ~typ:t.typ f' [trm_apps f [base']; index]
                | _ -> trm_map (aux global_trm) t
                end
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

let aos_to_soa (clog : out_channel) (name : var -> var) (x : typvar)
  (t : trm) : trm =
  (*
    changes:
    - declaration of x
    - add a copy of each function taking an argument of type x and replace the
    function calls with these copies
    - accesses on variables of type x
   *)
  let ilsm = functions_with_arg_type x t in
  (* first add copies of the functions *)
  let t = insert_fun_copies name ilsm x t in
  (* then replace function calls *)
  let t = replace_fun_names name ilsm x t in
  (* finally adapt the declaration and accesses *)
  swap_accesses clog x t


