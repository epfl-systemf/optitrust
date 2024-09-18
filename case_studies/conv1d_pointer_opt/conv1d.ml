open Optitrust
open Prelude
open Option

let _ = Flags.pretty_matrix_notation := false
let _ = Flags.disable_resource_typing ()

let hoist_if (if_tg: target) (tg : target) : unit =
  let pop_at (index : int) (t : trm) : trm * trm  =
    let error = "Sequence_core.delete_aux: expected the sequence on which the trms are deleted." in
    let tl = trm_inv ~error trm_seq_inv t in
    ((Mlist.nth tl index), trm_seq ~annot:t.annot (Mlist.remove index 1 tl)) in
  let inst: trm option ref = ref None in
  Target.iter (fun p ->
    let p, i = Internal.isolate_last_dir_in_seq p in
    apply_at_path (fun t ->
      let inst', seq = (pop_at i t) in
      let _ = inst := Some inst' in
      seq) p
  ) tg;
  Target.iter (fun p ->
    let p, i = Internal.isolate_last_dir_in_seq p in
    apply_at_path (fun t -> Sequence_core.insert_at (Option.get !inst) i t) p
  ) if_tg

let int_lit_from_trm (t: trm): int =
  let t_simpl = Arith_core.simplify false Arith_core.compute t in
  match t_simpl.desc with
    | Trm_lit (Lit_int i) -> i
    | _ -> trm_fail t "expected all constants in trm"

let inline_array_access (array_var : var) (new_vars: var list) (new_size : int) (t : trm) : trm =
  let rec aux (t : trm) : trm =
    match t.desc with
      | Trm_apps ({desc = Trm_prim (Prim_binop Binop_array_access);_} as p2, [{desc = Trm_var x;_}; index2], _) when (var_eq x array_var) ->
        let index_simpl = int_lit_from_trm index2 in
        let new_var = List.nth new_vars (index_simpl / new_size) in
        let new_index = trm_lit (Lit_int (index_simpl mod new_size)) in
        (trm_apps p2 [(trm_var new_var); new_index])
      | _ ->  trm_map aux t
    in aux t

let analyze_malloc (args:trm list) : int =
  match args with
    | _::sizes -> List.fold_left ( * ) 1 (List.map int_lit_from_trm (List.take ((List.length sizes) - 1) sizes))
    | _ -> failwith "not enough arguments to malloc"

let unroll_buffer_at (index : int) (t : trm) : trm =
  match t.desc with
  | Trm_seq tl ->
    let array_var = ref dummy_var in
    let new_size = ref 0 in
    let new_dims = ref 0 in
    let new_vars = ref [] in
    let process_decl (t : trm) : trm list =
      begin match t.desc with
        | Trm_let ((x , tx), init) ->
          array_var := x;
          let rec aux (t : trm): trm =
            begin match trm_apps_inv t with
            | Some ({desc = Trm_var f; _}, args) when (String.starts_with ~prefix:"MALLOC" f.name) ->
              let num_vars = int_lit_from_trm (List.hd args) in
              new_size := analyze_malloc args;
              array_var := x;
              new_vars := List.map new_var (List.init num_vars (fun i -> x.name ^ (string_of_int i)));
              new_dims := (List.length args) - 2;
              (trm_apps  ?loc:t.loc (trm_var (find_var_in_current_ast ("MALLOC" ^ (string_of_int !new_dims)))) (List.take_last (!new_dims + 1) args));
            | _ -> trm_map aux t;
            end in
            let init' = aux init in
            begin match !new_vars with
            | [] -> trm_fail t "Expected target to point to MALLOC declaration"
            | _ -> List.map (fun v -> (trm_let (v, tx) init')) !new_vars
            end
        | _ -> trm_fail t "Arrays_core.to_variables_aux: expected a variable declaration"
        end
      in

    let f_update_further (t : trm) : trm =
      inline_array_access !array_var !new_vars !new_size t
      in
    let new_stmts = process_decl (Mlist.nth tl index) in
    let new_tl = Mlist.update_at_index_and_fix_beyond ~delete:true index (fun t -> t) f_update_further tl in
    let new_tl' = Mlist.insert_sublist_at index new_stmts new_tl in
    let _, mfree_ind, mfree_args = (Mlist.nth (Mlist.filter (fun (a,b,c) -> a) (Mlist.mapi (fun i t -> begin match trm_apps_inv t with
          | Some ({desc = Trm_var f; _}, args) when (String.starts_with ~prefix:"MFREE" f.name) -> (true, i, List.tl args)
          | _ -> (false, i, [])
          end
      ) new_tl')) 0) in
    let mfree_stmts = List.map (fun v -> (trm_apps (trm_var (find_var_in_current_ast ("MFREE" ^ (string_of_int !new_dims)))) ((List.take !new_dims mfree_args) @ [trm_var v]))) !new_vars in
    let new_tl'' = Mlist.insert_sublist_at (mfree_ind) mfree_stmts (Mlist.remove (mfree_ind) 1 new_tl') in
    trm_seq ~annot:t.annot ?loc:t.loc new_tl''
  | _ -> t


(*
For a buffer with D dimensions NDxND-1x..N1
Unroll top dimension (outermost; D-1th) dimension of a buffer,
yielding N new variables of size NxD-1x..N1
Array accesses to buffer are scanned and replaced with the appropriate variable invocation; Mindex not supported for now.
If the index is not a literal constant, raise error

Target should be pointing to the declaration of the buffer
*)
let unroll_buffer (tg: target): unit =
  apply_at_target_paths_in_seq (unroll_buffer_at) tg

let elim_seq_surrounding (tg:target): unit =
  !! Marks_basic.add "seq_cln" tg;
  !! Sequence_basic.elim [cSeq ~args_pred:(Target.target_list_one_st [cMark "seq_cln"]) ()];
  !! Marks_basic.remove "seq_cln" tg


let nuke_nobraces (tg: target): unit =
  Target.iter (fun p ->
    let is_nobraces (t: trm): bool = List.fold_left (fun acc style -> begin match style with
        | No_braces _ -> acc || true
        | _ -> acc
      end) false t.annot.trm_annot_cstyle in
    let rec aux (t: trm): trm =
      begin match t.desc with
      | Trm_seq tl ->
        let ind = ref [] in
        let inst_list = ref [] in
        let tl' = Mlist.map (trm_map aux) tl in
        Mlist.iteri (fun i t ->
          let t' = trm_map aux t in
          if (is_nobraces t') then begin
            ind := i :: !ind;
            inst_list := (trm_inv trm_seq_inv t') :: !inst_list;
          end;
        ) tl;
        let _,tl'' = List.fold_left (fun (sum,tl_acc) (i,il) -> (sum+(Mlist.length il)-1), (Mlist.insert_sublist_at (sum+i) (Mlist.to_list il) (Mlist.remove (sum+i) 1 tl_acc))) (0,tl') (List.combine !ind !inst_list) in

        (trm_seq tl'')
      | _ -> trm_map aux t end
    in
    apply_at_path aux p
) tg

let extract_coeff (v: var) (t: trm): trm option =
  match (trm_binop_inv Binop_mul t) with
  | Some (lhs, rhs) ->
    let match_expr (lhs: trm) (rhs: trm): trm option = begin match (trm_var_inv rhs) with
      | Some(u) when (u.name = v.name) -> Some(lhs)
      | _ -> begin match (trm_var_get_inv rhs) with
          | Some(u) when (u.name = v.name) -> Some(lhs)
          | _ -> None
          end
        end in
      let l_match = match_expr lhs rhs in
      let r_match = match_expr rhs lhs in
      begin match l_match with
      | Some(t) -> l_match
      | None -> r_match
      end
  |_ -> begin match (trm_var_inv t) with
  | Some (u) when (u.name = v.name) -> Some (trm_lit (Lit_int 1))
  | _ -> None
  end

let insert_at_for_loop_end (vt: trm) (ft: trm): trm =
  begin match ft.desc with
  | Trm_for (lr,ts,lc) -> (trm_for ~annot:ft.annot ?loc:ft.loc ~contract:lc lr (trm_seq (Mlist.of_list [ts; vt])))
  | _ -> ft
  end

let trm_get (v:var): trm =
  (trm_apps (trm_unop Unop_get) [trm_var v])

let loopize_decl ?(loop_ind_in = dummy_var) (loop_tg: target) (decl_tg: target): unit =
  let loop_ind = ref loop_ind_in in
  if loop_ind_in = dummy_var then
    Target.iter_at_target_paths (fun t ->
      let l_range,_,_ = trm_inv trm_for_inv t in
      loop_ind := l_range.index;)
      loop_tg;
  let coeff = ref trm_dummy in
  let updated_decl = ref trm_dummy in
  Target.apply_at_target_paths (fun t ->
    let rec aux (t: trm): trm =
      let t' = match trm_apps_inv t with
      | Some ({desc = Trm_prim (Prim_unop Unop_get); _}, [a]) -> a
      | _ -> t in
      match extract_coeff !loop_ind t' with
      | Some c -> coeff := c; if (loop_ind_in = dummy_var) then trm_lit (Lit_int 0) else (trm_get loop_ind_in)
      | None -> trm_map aux t in
    updated_decl := aux t;
    if !coeff = trm_dummy then trm_fail t "can't extract a term with loop variable";
    trm_dummy
  ) decl_tg;
  apply_at_target_paths_in_seq (fun i t ->
    begin match t.desc with
    | Trm_seq tl ->
      let d_var,_,_ = trm_inv trm_let_inv !updated_decl in
      let acc_stmt = trm_apps (trm_binop Binop_set) [(trm_var d_var); (trm_add (trm_get d_var) !coeff)] in
      let tl' = Mlist.update_nth i (insert_at_for_loop_end acc_stmt) tl in
      let tl'' = Mlist.insert_at i !updated_decl tl' in
      trm_seq ~annot:t.annot ?loc:t.loc tl''
    | _ -> t
    end
  ) loop_tg

let loopize_expr ?(loop_ind_in = dummy_var) (loop_tg: target) (new_acc_name: string) (expr_tg: target): unit =
    let loop_ind = ref loop_ind_in in
    let new_acc_var = new_var new_acc_name in
    if loop_ind_in = dummy_var then
      Target.iter_at_target_paths (fun t ->
        let l_range,_,_ = trm_inv trm_for_inv t in
        loop_ind := l_range.index;)
        loop_tg;
    let coeff = ref trm_dummy in
    Target.apply_at_target_paths (fun t ->
      let rec aux (t: trm): trm =
        let t' = match trm_apps_inv t with
        | Some ({desc = Trm_prim (Prim_unop Unop_get); _}, [a]) -> a
        | _ -> t in
        match extract_coeff !loop_ind t' with
        | Some c -> coeff := c; trm_get new_acc_var
        | None -> trm_map aux t in
      let t' = aux t in
      if !coeff = trm_dummy then trm_fail t "can't extract a term with loop variable";
      t'
    ) expr_tg;
    apply_at_target_paths_in_seq (fun i t ->
      begin match t.desc with
      | Trm_seq tl ->
        let new_decl = (trm_let (new_acc_var,typ_ptr (typ_int)) (trm_apps (trm_prim (Prim_ref typ_int)) [trm_lit (Lit_int 0)])) in
        let acc_stmt = trm_apps (trm_binop Binop_set) [(trm_var new_acc_var); (trm_add (trm_get new_acc_var) !coeff)] in
        let tl' = Mlist.update_nth i (insert_at_for_loop_end acc_stmt) tl in
        let tl'' = Mlist.insert_at i new_decl tl' in
        trm_seq ~annot:t.annot ?loc:t.loc tl''
      | _ -> t
      end
    ) loop_tg

let int = trm_int

let cMindexAcc (arr_name: string): constr = (cBinop ~lhs:[cVar arr_name] ~rhs:[cFun ~regexp:true "MINDEX.*"] Binop_array_access)

let _ = Run.script_cpp (fun () ->
  !! Variable.bind_multi ~dest:[cIf (); tBefore] "ofs" [sExpr "j + r + 4 * tile_j"];
  !! nuke_nobraces [cFunDef "conv1d"];
  !! loopize_decl [cFor "r"; occFirst] [cVarDef "ofs"];
  !! Variable.bind "drow_ofs" [sExpr "c * 16 + ofs"];
  !! hoist_if [cIf ()] [cVarDef "drow_ofs"];
  !! loopize_decl ~loop_ind_in:(find_var_in_current_ast "ofs") [cFor "r"; occFirst] [cVarDef "drow_ofs"];
  !! loopize_expr [cFor "tile_j"] "data_base" [sExpr "4 * tile_j"];
  !! loopize_expr [cFor "c"] "data_row" [ sExpr "c * 16"]
)
