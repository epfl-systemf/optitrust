open Ast
open Target
open Ast_to_c
open Path_constructors
open Transformations
open Tools

(* loop_swap: Swap the two loop constructs, the loop should contain as least one innet loop
    params:
      path_to_loop an explicit path toward the loop
    returns: 
      the modified ast
*)

let rec loop_swap_core (clog : out_channel) (path_to_loop : path) (t : trm) =
  let (t,_) = resolve_path path_to_loop t in
  let log : string =
    let loc : string =
      match t.loc with 
      | None -> ""
      | Some (_,start_row,end_row,start_column,end_column) -> Printf.sprintf  "at start_location %d  %d end location %d %d" start_row start_column end_row end_column
    in 
    Printf.sprintf
    ("   - expression\n%s\n" ^^
    " %s  is a loop"
    )(ast_to_string t) loc 
    in write_log clog log;
    | Trm_labelled (l, t_loop) ->
    trm_labelled l (loop_swap_core clog t_loop)
    | Trm_seq [t_loop; t_del] when t.annot = Some Delete_instructions ->
      let t_swaped = loop_swap_core clog t_loop in
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
          | _ -> fail body1.loc "loop_swap_core: expected inner loop"
          end
      | _ -> fail t_swaped.loc "loop_swap_core: expected outer loop"
      end
    (* otherwise, just swap  *)
    | Trm_for (init1, cond1, step1,body1) ->
      let log : string =
        let loc : string =
          match body1.loc with
          | None -> ""
          | Some (_,start_row,end_row,start_column,end_column) -> Printf.sprintf  "at start_location %d  %d end location %d %d" start_row start_column end_row end_column
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
              | Some (_,start_row,end_row,start_column,end_column) -> Printf.sprintf  "at start_location %d  %d end location %d %d" start_row start_column end_row end_column
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
          | _ -> fail t.loc "loop_swap_core: inner_loop was not matched"
          end
        | _ -> fail t.loc "loop_swap_core; expected inner loop"
        end


    | _ -> fail t.loc "loop_swap_core; bad loop body"

    
(* loop_color: Replace the original loop with two nested loops:
        for (int i_color = 0; i_color < C; i_color++)
          for ([int] i = a+i_color; i < b; i += C) 
      params:
        path_to_loop: explicit path toward the loop
        c: a variable used to represent the number of colors
        i_color: index used for the new outer loop 
      returns:
        the modified ast
*)

let rec loop_coloring_core (clog : out_channel) (c : var) (new_var : var)(t : trm) : trm =
    let log : string =
      let loc : string =
        match t.loc with
        | None -> ""
        | Some (_,start_row,end_row,start_column,end_column) -> Printf.sprintf  "at start_location %d  %d end location %d %d" start_row start_column end_row end_column
      in
      Printf.sprintf
        ("  - expression\n%s\n" ^^

         "    %sis a (labelled) loop\n"
        )
        (ast_to_string t) loc
      in
    write_log clog log;
    match t.desc with
    (* The loop might be labelled, so keep the label *)
    | Trm_labelled (l, t_loop) ->
      trm_labelled l (loop_coloring_core clog c new_var t_loop)

    | Trm_seq [t_loop; t_del] when t.annot = Some Delete_instructions ->
     let t_transformed = loop_coloring_core clog c new_var t_loop in
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
        | _ -> fail body1.loc "loop_coloring_core: expected inner loop"
        end
     | _ -> fail t_transformed.loc "loop_coloring_core: expected outer loop"
     end
  | Trm_for (init, cond, step, body) ->
      let log : string =
        let loc : string =
          match body.loc with
          | None -> ""
          | Some (_,start_row,end_row,start_column,end_column) -> Printf.sprintf  "at start_location %d  %d end location %d %d" start_row start_column end_row end_column
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


  | _ -> fail t.loc "loop_coloring_core: not a for loop, check the path "

