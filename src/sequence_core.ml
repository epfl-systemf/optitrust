open Ast
(* *********************************************************************************** 
 * Note: All the intermediate functions which are called from [sequence.ml] file      *
 * have only one purpose, and that is targeting the trm in which we want to apply the *
 * transformation. That's why there is not need to document them.                     *
 *)

(* [insert_aux index ts t]: insert an arbitrary trm after or before a targeted trm
    params:
      index: and integer in range 0 .. (nbinstr-1)
      ts: a list of trms which the inner sequence will contain
      t: ast of the outer sequence where the insertion will be performed.
    return: 
      the updated ast of the updated outer sequence with the augmented new trm
*)
let insert_aux (index : int) (s : string) (t : trm) : trm =
    match t.desc with
    | Trm_seq tl ->
      let new_trm = code s in
      let new_tl = Mlist.insert_at index new_trm tl in
      trm_seq ~annot:t.annot ~marks:t.marks new_tl
    | _ -> fail t.loc "insert_aux: expected the sequence on which the insertion is performed"

let insert (index : int) (s : string) : Target.Transfo.local =
  Target.apply_on_path (insert_aux index s)

(* [delete_aux index nb_instr t]: delete a number of instruction inside the sequence starting 
      from index [index] and ending at ([index] + [nb])
    params:
      nb: number of instructions to delete
      t: ast of the outer sequence where the deletion is performed.
    return: 
      the updated ast of the surrounding sequence with the deleted nodes
*)
let delete_aux (index : int) (nb_instr : int) (t : trm) : trm =
  match t.desc with
    | Trm_seq tl ->
      trm_seq ~annot:t.annot ~marks:t.marks (Mlist.remove index nb_instr tl)
    | _ -> fail t.loc "delete_aux: expected the sequence on which the trms are deleted"

(* [delete index nb_instr t p] *)
let delete (index : int) (nb_instr : int) : Target.Transfo.local =
  Target.apply_on_path (delete_aux index nb_instr)


(* [intro_aux index nb t]: inside a sequence, move all the trms with findex falling in a rance 
      from [index] to [index] + [nb] into a sub-sequence.
    params:
      index: index where the grouping is performed
      ts: a list of ast nodes
      t: ast of the outer sequence where the insertion is performed
    return: the updated ast of surrounding sequence with the inserted nodes

*)

let intro_aux (mark : string) (index : int) (nb : int) (t : trm) : trm =
  match t.desc with
    | Trm_seq tl ->
      let tl1, tl2 = 
        if nb > 0 then Mlist.extract index nb tl else Mlist.extract (index+ nb+1) (-nb) tl in
        let intro_seq = trm_seq tl2 in
        let intro_seq = if mark <> "" then trm_add_mark mark intro_seq else intro_seq in
        let index = if nb < 0 then index -1 else index in
         trm_seq  ~annot:t.annot ~marks:t.marks (Mlist.insert_at index intro_seq tl1)
    | _ -> fail t.loc "intro_aux: expected the sequence on which the grouping is performed"

let intro (mark : string) (index : int) (nb_instr : int) : Target.Transfo.local =
  Target.apply_on_path (intro_aux mark index nb_instr)

(*[elim_aux index t]: inline an inner sequence into an outer sequence.
    params:
      t: ast of the sequence wanted to remove
    return: 
      a hiden sequence which is going to be merged witht the outer sequence on the next step
*)
let elim_aux (t : trm) : trm =
  match t.desc with
  | Trm_seq tl ->
     trm_seq_no_brace (Mlist.to_list tl)
  | _ -> fail t.loc "elim_aux: expected the sequence to be deleteds"

let elim : Target.Transfo.local =
  Target.apply_on_path(Internal.apply_on_path_targeting_a_sequence (elim_aux) "elim")

(* [intro_on_instr_aux visible mark t]: replacing t with a sequence that contains t as single item.
   params:
    mark: add a mark around the sequence
    visible: a flag to turn on(off) curly braces of the sequence
    t: ast of the instruction 
   return: 
    updated ast of the outer sequence with wrapped node t
 *)
let intro_on_instr_aux (mark : mark) (visible : bool) (t : trm) : trm =
  let wrapped_seq = if visible then trm_seq (Mlist.of_list [t]) else trm_seq_no_brace [t] in
  if mark <> "" then trm_add_mark mark wrapped_seq else wrapped_seq 
 
let intro_on_instr (visible : bool) (mark : mark) : Target.Transfo.local=
  Target.apply_on_path (intro_on_instr_aux mark visible)

(* [unrwap_aux t]: replacing a sequence that contains a single item t with t.
   params:
    t: a term that corresponds to a sequence with a single item in t
   return:
    the udated the ast where the trm inside the sequence has been extracted
 *)
let unwrap_aux (t : trm) : trm =
  match t.desc with
    | Trm_seq tl ->
      if Mlist.length tl = 1 then Mlist.nth tl 0 
        else fail t.loc "unwrap_aux: can only unwrap a sequence with exactly one item"
    | _ -> fail t.loc "unwrap_aux: expected to operate on a sequence"

let unwrap : Target.Transfo.local =
  Target.apply_on_path (unwrap_aux)

(* [split_aux index t ]: splitting a sequence in two parts 
    params:
      index: index of the realative target entered by the user
      t : an ast term that corresponds to the the targeted sequence
    return:
      the updated ast with the splitted sequence
*)
let split_aux (index : int) (t : trm) : trm =
  match t.desc with 
  | Trm_seq tl ->
    let first_part,last_part = Mlist.split index tl in
    trm_seq_no_brace [trm_seq ~annot:t.annot first_part;trm_seq ~annot:t.annot last_part]
  | _ -> fail t.loc "split_aux: expected a sequence, containing the location where it is going to be splitted"

let split (index : int) : Target.Transfo.local =
  Target.apply_on_path (split_aux index)


let partition_aux (blocks : int list) (braces : bool) (t : trm) : trm =
  Tools.printf "%b\n" (Internal.is_nobrace t);
  match t.desc with 
  | Trm_seq tl -> 
    let nb = Mlist.length tl in
    let blocks = if blocks = [] then [nb] else blocks in
    let sum_blocks = List.fold_left (+) 0 blocks in
    if sum_blocks <> nb 
      then fail t.loc (Tools.sprintf "partition: the partition entered is not correct, the list length is %d, while the sum of the block size is %d" (Mlist.length tl) sum_blocks)
      else
        let current_list = ref tl in
        let partition = List.fold_left (fun acc x -> 
          let lfront, lback = Mlist.split x !current_list in
          current_list := lback;
          lfront :: acc
          ) [] blocks 
          in
        let new_tl = 
          if braces 
            then Mlist.of_list (List.map trm_seq (List.rev partition))
            else Mlist.of_list (List.map (fun x -> trm_seq_no_brace (Mlist.to_list x)) (List.rev partition))
            in
        trm_seq ~annot:t.annot ~marks:t.marks new_tl
        
  | _ -> fail t.loc "partial_aux: expected a sequence to partition"

let partition (blocks : int list) (braces : bool): Target.Transfo.local =
  Target.apply_on_path (partition_aux blocks braces)


let shuffle_aux (braces : bool) (t : trm) : trm =
  match t.desc with 
  | Trm_seq tl ->
    if Mlist.length tl < 1 then fail t.loc "shuffle_aux:can't shuffle an empty mlist";
    let first_row = Mlist.nth tl 0 in
    begin match first_row.desc with 
    | Trm_seq tl1 ->
      let loop_bound = Mlist.length tl1 in
      if loop_bound < 2 then fail t.loc "shuffle_aux: expected a row of length at least 2";
      let global_acc = ref [] in
      for i = 0 to loop_bound-1 do
        let local_acc = Mlist.fold_left (fun acc t1 -> 
            begin match t1.desc with 
            | Trm_seq tl2 ->
              if Mlist.length tl2 <> loop_bound then fail t1.loc "shuffle_aux: all the subgroups should be of the same size";
              let temp_el = Mlist.nth tl2 i in
              let temp_el = 
              if braces 
                then Internal.remove_nobrace_if_sequence temp_el 
                else Internal.set_nobrace_if_sequence temp_el in
            temp_el :: acc
            | _ -> fail t1.loc "shuffle_aux: all the elements of the blocks should be sequences"
            end
            
          ) [] tl in
        let local_acc = List.rev local_acc in
        global_acc := (if braces then trm_seq (Mlist.of_list local_acc) else trm_seq_no_brace local_acc) :: !global_acc
      done;
       trm_seq ~annot:t.annot ~marks:t.marks (Mlist.of_list (List.rev !global_acc))

    | _ -> fail first_row.loc "shuffle_aux: shuffle can be applied only on sequences"
    end
  | _ -> fail t.loc "shuffle_aux: expected the sequence with blocks to reorder"

let shuffle (braces : bool) : Target.Transfo.local = 
  Target.apply_on_path (shuffle_aux braces) 