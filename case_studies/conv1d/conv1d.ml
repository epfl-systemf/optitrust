open Optitrust
open Prelude
open Option

let _ = Flags.pretty_matrix_notation := false
let _ = Flags.disable_resource_typing ()

let replace_var_name (new_var_name: string) (tg : target) : unit =
  let replace_var_name_trm (new_var_name : string) (t: trm) : trm =
    let error = "Expected a target to a variable" in
    let v = trm_inv ~error trm_var_inv t in
    (* to change name and keep namespaces/id:
    let fv = trm_inv ~error trm_var_inv f in
    { namespaces = fv.namespaces; name = new_fun_name; id = fv.id } *)
    trm_replace (Trm_var {id = v.id; name = new_var_name; namespaces = v.namespaces}) t in
    apply_at_target_paths (replace_var_name_trm new_var_name) tg

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

let if_split (if_tg: target) (tg: target) : unit =
  let split_at (index : int) (t : trm) : trm * trm  =
    let error = "Sequence_core.delete_aux: expected the sequence on which the trms are deleted." in
    let tl = trm_inv ~error trm_seq_inv t in
    let lhs, rhs = (Mlist.split (index+1) tl) in
    (trm_seq ~annot:t.annot lhs, trm_seq ~annot:t.annot rhs) in
  let rhs: trm option ref = ref None in
  Target.iter (fun p ->
    let p, i = Internal.isolate_last_dir_in_seq p in
    apply_at_path (fun t ->
      let seq, rhs' = (split_at i t) in
      let _ = rhs := Some rhs' in
      seq) p
  ) tg;
  let cond,_,_ = trm_inv trm_if_inv (Option.get (get_trm_at if_tg)) in
  Target.iter (fun p ->
    let p, i = Internal.isolate_last_dir_in_seq p in
    apply_at_path (fun t -> Sequence_core.insert_at (trm_if (cond) (Option.get !rhs) (trm_unit ())) (i+1) t) p
  ) if_tg

let debug_loop (tg : target) : unit =
  Target.iter (fun p ->
    let tg_loop = Target.resolve_path p in
    let indices = Internal.get_loop_nest_indices tg_loop in
    print_string (String.concat " " (List.map (fun v -> v.name) indices)))
    tg

(* TODO assert that number of dims has to be >= 2*)
let reduce_mindex (dim: int) (app_target: target): unit =
    Target.iter (fun p ->
        apply_at_path (fun t ->
          let _,args = trm_inv trm_apps_inv t in
          let dims = List.length args / 2 in
          let lengths, indices = List.split_at dims args in
          let hi_term = (List.nth lengths (dims - (dim+2))) in
          let hi_term_m = if (dim+2 = dims) then
            trm_add_mark "window_removable" hi_term
          else hi_term in
          let lengths' = (List.take (dims - (dim+2)) lengths) @ [trm_mul (List.nth lengths (dims - (dim+1))) hi_term_m] @ (List.take_last dim lengths) in
          let indices' = (List.take (dims - (dim+2)) indices)
           @ [trm_add (trm_mul (List.nth lengths (dims - (dim+1))) (List.nth indices (dims - (dim+2)))) (List.nth indices (dims - (dim+1)))]
           @ (List.take_last dim indices) in
          trm_apps (trm_var (find_var_in_current_ast ("MINDEX"^(string_of_int (dims-1))))) (lengths' @ indices')
        ) p
    ) app_target

let is_prim_add (p : prim) : bool =
  match p with
  | Prim_binop (Binop_add) -> true
  | _ -> false

let extract_window_offset (mindex_trm: trm) (keep_indices: string list): trm * trm =
  let f,args = trm_inv trm_apps_inv mindex_trm in
  let dims = List.length args / 2 in
  let lengths, indices = List.split_at dims args in
  let indices_enum = List.mapi (fun i a -> (i,a)) indices in
  let adjust_strides (a: trm) (ind: int): trm =
    let strides = List.take_last (dims - ind - 1) lengths in
    match strides with
      | [] -> a
      | _ -> List.fold_left (fun acc s -> trm_mul acc s) a strides in
  let offsets, args = List.fold_left (fun (acc: trm list * trm list) (ind,a) -> match (trm_add_inv a) with
   | Some (lhs, ({desc=Trm_var v;_} as rhs)) when (List.mem v.name keep_indices) -> (((adjust_strides lhs ind) :: (fst acc)), (snd acc) @ [rhs])
   | Some (({desc=Trm_var v;_} as lhs), rhs) when (List.mem v.name keep_indices) -> (((adjust_strides rhs ind) :: (fst acc)), (snd acc) @ [lhs])
   | Some (_,_) -> ((adjust_strides a ind) :: (fst acc), snd acc)
   | None -> begin match (trm_var_inv a) with
    | Some (v) when (List.mem v.name keep_indices) -> (fst acc, (snd acc) @ [a])
    | _ -> (((adjust_strides a ind) :: (fst acc)), (snd acc))
   end
  ) ([],[]) indices_enum in
  let lengths' = List.map (fun l -> match (trm_binop_inv Binop_mul l) with
   | Some (a,b) when (trm_has_mark "window_removable" a) -> b
   | Some (a,b) when (trm_has_mark "window_removable" b) -> a
   | _ -> l
  ) lengths in
  let offset = List.fold_left (fun acc a -> (trm_add acc a)) (List.hd offsets) (List.tl offsets) in
  (offset,(trm_apps f (lengths' @ args)))


let factor_window (keep_indices: string list) (access_target: target): unit =
  Target.iter (fun p ->
    apply_at_path (fun t ->
      let f,args = trm_inv ~error:"not an array access?" trm_apps_inv t in
      (* TODO assert the prim is array access *)
      let arr = List.nth args 0 in
      let access_expr = List.nth args 1 in
      let offset, rem_access_expr = extract_window_offset access_expr keep_indices in
      trm_apps (trm_prim (Prim_binop Binop_array_access)) [(trm_apps (trm_prim (Prim_binop Binop_array_access)) [arr; offset]); rem_access_expr]
      (*trm_apps (trm_prim (Prim_binop Binop_array_access)) [(trm_add arr offset); rem_access_expr]*)

    ) p
) access_target

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

(*
 begin match (trm_var_inv rhs) with
      | Some(u) when (u.name = v.name) -> Some(c)
      | _ -> begin match (trm_var_get_inv rhs) with
          | Some(u) when (u.name = v.name) -> Some(c)
          | _ -> None
          end
        end*)
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
  !! Loop.tile (trm_int 4) ~index:"tile_i" ~iter:TileIterGlobal ~bound:TileDivides [cFor "i"];
  !! Loop.tile (trm_int 4) ~index:"tile_j" ~iter:TileIterGlobal ~bound:TileDivides [cFor "j"];
  !! Loop.reorder ~order:["tile_j"; "i"] [cFor "i"];
  !! Loop.tile (trm_int 4) ~index:"tile_i_hi" ~iter:TileIterGlobal ~bound:TileDivides [cFor "tile_i"];
  !! Loop.reorder ~order:["tile_j"; "tile_i"] [cFor "tile_i"];

  !! Loop.hoist_alloc_loop_list [1; 1; 1;] [cVarDef "sum"];
  (* zeroing *)
  !! Loop.fission ~nest_of:3 [tBefore; cFor "k"; occFirst];
  (* storing *)
  !! Loop.fission ~nest_of:3 [tAfter; cFor "k"; occFirst];

  let prefix_indices (idxs: (string * target) list) (pfx: string): unit =
    let _ = (List.map (fun idx -> Loop.rename_index (pfx ^ "_" ^ fst idx) (snd idx)) idxs) in () in
  !! prefix_indices [("tile_i", [cFor "tile_i"; occFirst]); ("i", [cFor "i"; occFirst]); ("j", [cFor "j"; occFirst])] "zero";
  !! prefix_indices [("tile_i", [cFor "tile_i"; occLast]); ("i", [cFor "i"; occLast]); ("j", [cFor "j"; occLast])] "st";

  !! Loop.reorder ~order:["k"; "j"] [cFor "j"];
  !! Loop.reorder ~order:["k"; "i"] [cFor "i"];
  !! Loop.reorder ~order:["k"; "tile_i"] [cFor "tile_i"];

  (* TODO: interesting that I can get these wrong and the program still compiles? *)
  (* The expression should be dependent on j and r, but i can make it independent of r*)
  !! Loop.hoist_alloc_loop_list [0; 0; 1; 1;] [cVarDef "y"];
  !! if_split [cIf ();] [cArrayWrite "y"; occLast];
  !! Loop.fission ~nest_of:4 [tAfter; cIf (); occFirst];


  let delete_loop (loop_tg: target) (inner_loop_tg: target): unit =
    Loop_basic.move_out inner_loop_tg;
    Loop.delete_void loop_tg in

  delete_loop [cFor "tile_i"; occFirst] [cFor "i"; occFirst];
  delete_loop [cFor "i"; occFirst] [cFor "j"; occFirst];

  !!hoist_if [cIf (); occLast] [cArrayWrite "sum"; occLast];
  !!Instr_basic.delete [cIf(); occLast];

  (* BEGIN RVM OPTIMIZATIONS *)
  (*let stage_mem ?(extra_transforms: string -> unit = fun _ -> ()) (name: string) (nest_levels: int) (expr_target: target_simple): unit =
  Variable.bind name expr_target;
  extra_transforms name;
  Loop.fission ~nest_of:nest_levels [tAfter; cVarDef name];*)

  !! Variable.bind "data_tile" [cArrayRead "y"];
  Loop.hoist_alloc_loop_list [0; 0; 1; 1;] [cVarDef "data_tile"];
  Loop.fission ~nest_of:4 [tAfter; cArrayWrite "data_tile"];

  !! delete_loop [cFor "tile_i"; occFirst] [cFor "i"; occFirst];
  delete_loop [cFor "i"; occFirst] [cFor "j"; occIndex 1];

  !! Function.uninline ~contains_for_loop:false ~fct:[cFunDef "rvm_mld"] [cFor "j"; occIndex 1];

  !! Variable.bind "kernel_tile" [cArrayRead "W"];
  Loop.hoist_alloc_loop_list [1; 0; 1;] [cVarDef "kernel_tile"];
  Loop.fission ~nest_of:3 [tAfter; cArrayWrite "kernel_tile"];
  !! delete_loop [cFor "j"; occIndex 1] [cFor "r"; occIndex 1];

  !! reduce_mindex 0 [cFun "MINDEX3"; occIndex 1];
     factor_window ["i"; "r"] [cMindexAcc "W";];
  !!   Function.uninline ~contains_for_loop:false ~fct:[cFunDef "rvm_mld"] [cFor "i"; occIndex 0];

  !! reduce_mindex 1 [cFun "MINDEX3"; occIndex 1];
    factor_window ["i"; "j"] [cMindexAcc "sum"; occIndex 1];
  (* TODO: had to commute the multiply in the spec? surely there is a simple rewrite for this? *)
  !!  Function.uninline ~contains_for_loop:false ~fct:[cFunDef "rvm_mma"] [cFor "i"; occIndex 0];

  !! reduce_mindex 1 [cFun "MINDEX3"; occIndex 1];
  !!   factor_window ["st_i"; "st_j"] [cMindexAcc "sum"; occIndex 1];
  !! factor_window ["st_i"; "st_j"] [cMindexAcc "O";];
  !!  Function.uninline ~contains_for_loop:false ~fct:[cFunDef "rvm_mst"] [cFor "st_i"];

  !! reduce_mindex 1 [cFun "MINDEX3"; occIndex 0];
  !!   factor_window ["zero_i"; "zero_j"] [cMindexAcc "sum"; occIndex 0];
  !! Function.uninline ~contains_for_loop:false ~fct:[cFunDef "rvm_mzero"] [cFor "zero_i"];

  !! Loop.unroll [cFor "zero_tile_i"];
  !! Loop.unroll [cFor "st_tile_i"];
  !! Loop.unroll [cFor "tile_i"];
  !! unroll_buffer [cVarDef "sum"];

  !! Variable.bind_multi ~dest:[cIf (); tBefore] "ofs" [sExpr "tile_j * 4 + j + r"];
  !! nuke_nobraces [cFunDef "conv1d"];
  !! Matrix_basic.elim_mindex [cMindex ~args:[[cVar "IC"]; [cVar "IW"]; [cVar "k"]; [cVar "ofs"]] ()];
  !! Matrix_basic.elim_mindex [nbMulti; cMindex ~args:[[cInt 4]; [cInt 4]; [cVar "j"]; [cVar "r"]] ()];

  !! loopize_decl [cFor "r"; occFirst] [cVarDef "ofs"];
  !! Arith_basic.simplify [sExpr "0 + k * IW + ofs"];
  !! Variable.bind "drow_ofs" [sExpr "k * IW + ofs"];
  !! hoist_if [cIf ()] [cVarDef "drow_ofs"];
  !! loopize_decl ~loop_ind_in:(find_var_in_current_ast "ofs") [cFor "r"; occFirst] [cVarDef "drow_ofs"];
  !! loopize_expr [cFor "tile_j"] "data_base" [sExpr "tile_j * 4"];
  !! loopize_expr [cFor   "k"] "data_row" [sExpr "k * IW"]
  (*   Variable.bind_multi ~dest:[cArrayWrite "y"; tBefore; occFirst] "dtile_ofs" [sExpr "0 + j * 4 + r"];
   nuke_nobraces [cFunDef "conv1d"];*)
)

  (*!! Variable_basic.init_detach [cVarDef "y"];
  !! Instr_basic.copy  [cWriteVar "y"];
  !! Expr_basic.replace (stmt "y = 0") [cWriteVar "y"; occFirst];
  !! Variable_basic.init_attach [cVarDef "y"];*)
  (* TODO: this step is only _sort of_ valid; it's not going to cause incorrect behavior,
     unless the access out of bounds crashes the program *)
  (*!! hoist_if [cIf ()] [cVarDef "y"];*)
  (*~dest:[tBefore; cFor "i"]*)
