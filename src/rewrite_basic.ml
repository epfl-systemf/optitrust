open Ast

(* [parse_pattern str]: for a given pattern [str] return the list of variables used in that pattern 
      and the ast of the pattern
*)
let parse_pattern (str : string) : (vars * trm) = 
  let output_file = "tmp_rule.cpp" in
  let file_content = "int main() {" ^ str ^";}" in
  Xfile.put_contents output_file file_content;
  let _, ast_of_file = Trace.parse output_file in
  match ast_of_file.desc with 
  | Trm_seq tl when (List.mem Main_file ast_of_file.annot) -> 
    if Mlist.length tl = 0 then fail ast_of_file.loc "parse_pattern: couldn't parse pattern";
    let main_fun = Mlist.nth tl 0 in
    begin match main_fun.desc with 
    | Trm_let_fun (_, _, _, body) ->
      begin match body.desc with 
      | Trm_seq tl1 -> 
        let tl_vars, tl_pattern = Mlist.split (Mlist.length tl1 - 1) tl1 in                
        let pattern_vars = List.flatten (List.map trm_vardef_get_vars (Mlist.to_list tl_vars)) in
        let pattern_instr = Mlist.nth tl_pattern 0 in
        (pattern_vars, pattern_instr)
      | _ -> fail body.loc "parse_pattern: body of the main function should be a sequence"
      end
    | _ -> fail main_fun.loc "parse_pattern: expected the declaration of the main function"
    end
  | _ -> fail ast_of_file.loc "parse_pattern: expected the main sequence of tmp_rule.cpp"



(* [parse_rule str]: for a given pattern [str] return a rewrite rule which is a record containing the 
    the list of variables used in that rule, the rule itself and the result after applying that rule.
*)
let parse_rule (str : string) : rewrite_rule =
  let pattern_vars, pattern_instr = parse_pattern str in
  match pattern_instr.desc with 
  | Trm_apps ({desc = Trm_val (Val_prim (Prim_binop Binop_eq));_},[t1; t2]) ->
    {rule_vars = pattern_vars; rule_from = t1; rule_to = t2}
  | _ -> fail pattern_instr.loc "parse_rule: could not parse the given rule"

exception Rule_mismatch 

(* [rule_match]: for a given list of variables [vars] in pattern [pat] return a map with all the variables  
      as keys and their associated ast as values
*)
let rule_match (vars : vars) (pat : trm) (t : trm) : tmap =
  Tools.printf "Comparing %s with %s\n" (Ast_to_c.ast_to_string pat) (Ast_to_c.ast_to_string t);
  let inst = ref Trm_map.empty in
  let rec aux (t1 : trm) (t2 : trm) : unit =
    let aux_list (ts1 : trm list) (ts2 : trm list) : unit =
      List.iter2 aux ts1 ts2 in  
    match t1.desc, t2.desc with 
    (* | Trm_apps (_, [t']), _ when is_get_operation t1 -> aux t' t2
    | _, Trm_apps (_, [t']) when is_get_operation t2 -> aux t1 t' *)
    | Trm_var x, _ when List.mem x vars ->
      begin match Trm_map.find_opt x !inst with 
      | None -> inst := Trm_map.add x t2 !inst
      | Some t0 when (Internal.same_trm t0 t2) -> raise Rule_mismatch
      | _ -> ()
      end
    | Trm_var x1, Trm_var x2 when x1 = x2 -> ()
    | Trm_val v1, Trm_val v2 when Internal.same_val v1 v2 -> ()
    | Trm_apps (f1, ts1), Trm_apps (f2, ts2) ->
      begin 
        aux f1 f2;
        aux_list ts1 ts2;
      end
    | _ -> 
      raise Rule_mismatch
  in 
  aux pat t;
  !inst

(* [apply_rule rule t] apply rule [rule] in the ast [t] *)
let apply_rule (rule : rewrite_rule) (t : trm) : trm =
  let inst : tmap = rule_match rule.rule_vars rule.rule_from t in
  Internal.variable_substitute inst rule.rule_to


