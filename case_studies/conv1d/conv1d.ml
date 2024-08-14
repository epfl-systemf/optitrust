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
          let lengths' = (List.take (dims - (dim+2)) lengths) @ [trm_mul (List.nth lengths (dims - (dim+1))) (List.nth lengths (dims - (dim+2)))] @ (List.take_last dim lengths) in
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
  let offsets, args = List.fold_left (fun (acc: trm list * trm list) a -> match (trm_add_inv a) with
   | Some (lhs, ({desc=Trm_var v;_} as rhs)) when (List.mem v.name keep_indices) -> ((lhs :: (fst acc)), (snd acc) @ [rhs])
   | Some (({desc=Trm_var v;_} as lhs), rhs) when (List.mem v.name keep_indices) -> ((rhs :: (fst acc)), (snd acc) @ [lhs])
   | Some (_,_) -> (a :: (fst acc), snd acc)
   | None -> begin match (trm_var_inv a) with
    | Some (v) when (List.mem v.name keep_indices) -> (fst acc, (snd acc) @ [a])
    | _ -> ((a :: (fst acc)), (snd acc))
   end
  ) ([],[]) indices in
  let offset = List.fold_left (fun acc a -> (trm_add acc a)) (List.hd offsets) (List.tl offsets) in
  (offset,(trm_apps f (lengths @ args)))


let factor_window (keep_indices: string list) (access_target: target): unit =
  Target.iter (fun p ->
    apply_at_path (fun t ->
      let f,args = trm_inv ~error:"not an array access?" trm_apps_inv t in
      (* TODO assert the prim is array access *)
      let arr = List.nth args 0 in
      let access_expr = List.nth args 1 in
      let offset, rem_access_expr = extract_window_offset access_expr keep_indices in
      trm_apps (trm_prim (Prim_binop Binop_array_access)) [(trm_add arr offset); rem_access_expr]
    ) p
) access_target

let int = trm_int

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
  !! factor_window ["i"; "r"] [cBinop ~lhs:[cVar "W"] ~rhs:[cTrue] Binop_array_access];
  !! Function.uninline ~contains_for_loop:false ~fct:[cFunDef "rvm_mld"] [cFor "i"; occIndex 0];
)

  (*!! Variable_basic.init_detach [cVarDef "y"];
  !! Instr_basic.copy  [cWriteVar "y"];
  !! Expr_basic.replace (stmt "y = 0") [cWriteVar "y"; occFirst];
  !! Variable_basic.init_attach [cVarDef "y"];*)
  (* TODO: this step is only _sort of_ valid; it's not going to cause incorrect behavior,
     unless the access out of bounds crashes the program *)
  (*!! hoist_if [cIf ()] [cVarDef "y"];*)
  (*~dest:[tBefore; cFor "i"]*)
