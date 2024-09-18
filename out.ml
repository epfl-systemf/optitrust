
(* //////////////////////////////
   ////  Sequence_basic.ml  /////
   ////////////////////////////// *)

      (*========================*)
(** [insert ~reparse code tg]: expects the target [tg] to point at a relative position(in between two instructoins),
     [code] - the instruction that is going to be added, provided by the user as an arbitrary trm. *)
(*========================*)
(** [delete index nb tg]: expects the target [tg] to point at an instruction,
     [nb] - denotes the number of instructions to delete starting from the targeted trm.

   @correctness: correct if nothing modified by the instruction was observed
   later.
   If the next instructions need an invariant H' and { H } del_instr { H'' }
   we need both H ==> H' and H'' ==> H'. *)
(*========================*)
(** [intro i nb tg]: expects the target [tg] to point at an instruction inside a sequence.
    [mark] -  denotes a mark which add into the generated sub-sequence, in case the user decides to have one.
    [visible] - denotes the visibility of a sequence. This means the that the the sequence is
               used only for internal purposes.                     }
    [nb] - is the number of instructions to be moved inside the sub-sequence.
          If [nb] = 1 means then this transformation is basically the same as intro_on_instr.
          If [nb] is greater than one then it means that the instructions which come right after
            the targeted instruction will be included in the sub-sequence too.
          If [nb] is lwoer than one then it means that the instructions which come before
            the targeted instruction will be included in the sub-sequence too.
    Ex: int main(){     int main(){
        int x = 5;      { int x = 5}
        iny y = 6;      int y = 6;
        return 0;       return 0;
      }                } *)
(*========================*)
(** [intro_after ~mark ~label tg]: same as [intro] but this transformation will include in the sequence all the
    instructions that come after the targeted instruction and belong to the same scope. *)
(*========================*)
(** [intro_before ~mark ~label tg]: similar to [intro] but this transformation will include in the sequence all the
    instructions that come before the targeted instruction and belong to the same scope. *)
(*========================*)
(** [intro_between ~mark ~label tg_beg tg_end]: this transformation is an advanced version of [intro].
     Here, the user can specify explicitly the targets to the first and the last instructions that
     are going to be isolated into a sequence. *)
(*========================*)
(** [elim tg]: expects the target [tg] to point at a sequence that appears nested inside another sequence,
    e.g., points at [{t2;t3}] inside [{ t1; { t2; t3 }; t4 }]. It "elims" the contents of the inner sequence,
    producing e.g., [{ t1; t2; t3; t3}]. *)
(*========================*)
(** [intro_on_instr ~mark ~visible tg]: expecets the target [tg] to point at an instruction,
    then it will wrap a sequence around that instruction.
    [visible] - denotes the visibility of a sequence. This means the that the the sequence is
        used only for internal purposes.
    [mark] - denotes the mark of the sub-sequence. Targeting sequences can be challanging hence having
          them marked before can make the apllication of the transformations easier. *)
(*========================*)
(** [elim_on_instr tg]: expects the target [tg] to point at a sequence that contains a single instruction,
    then it removes that sequence. *)
(*========================*)
(** [split tg]: expects the target [tg] to point in between two instructions, then it will split the sequence
     that contains that location into two sequences. *)
(*========================*)
(** [partition ~braces blocks tg]: expects the target tg to point at a sequence, this transformations will split that sequence
      into blocks where the sizes of the blocks should be provided by the user.
        [blocks] - denotes the sizes for each block inside the sequence. By default it is empty, otherwise the sum of
          integers inside [blocks] should sum up to the number of instructions of the targeted sequence.
        [braces] - denotes a flag for the visibility of the blocks meaning that this block partition will be meaningful only for
          other transformations that call explicitly the partition transformation. *)
(*========================*)
(** [shuffle ~braces tg]: expects the target [tg] to point at a sequence of blocks, this transformation will transpose the block structure

    think about a sequence of blocks as a matrix.
    {
      {{t11};{t12};{t13}};
      {{t21};{t22};{t23}};
      {{t31};{t32};{t33}};
    }
    this will be changed to:
    {
      {{t11};{t21};{t31}};
      {{t12};{t22};{t32}};
      {{t13};{t23};{t33}};
    } *)

(* //////////////////////////////
   //////  Ghost_pure.ml  ///////
   ////////////////////////////// *)

      
(* //////////////////////////////
   //////  Marks_basic.ml  //////
   ////////////////////////////// *)

      (*========================*)
(** [add m tg]: adds mark [m] to the term or interstice that correpsonds to target [tg].
   NOTE: if m = "" then does nothing. *)
(*========================*)
(** [remove m tg]: removes mark m from the term or interstice that corresponds to target [tg]. *)
(*========================*)
(** [remove_st pred tg]: cleans all the marks satisfying [pred] from the trm
   that corresponds to target [tg]. Use [~indepth:true] to recurse in subterms. *)
(*========================*)
(** [clean tg]: removes all the marks form the term that is matched by target [tg].
   Use [~indepth:true] to recurse in subterms. *)

(* //////////////////////////////
   ///////  Omp_basic.ml  ///////
   ////////////////////////////// *)

      
(* //////////////////////////////
   ///////  Function.ml  ////////
   ////////////////////////////// *)

      (*========================*)
(** [rename]: instantiation of Rename module *)
(*========================*)
(** [bind_args fresh_names tg]: expects the target [tg] to point at a function call.
      Then it takes [fresh_names] which is a list of strings where the string
      at index i represents the variable going to be binded to the argument i
      of the function call. If one doesn't want to bind the argument at index i
      then it just leaves it as an empty string "". Basically this transformation is
      just an aplication of bind_intro n times. Where n is the numer of strings inside
      [fresh_names] different from "". *)
(*========================*)
(** [elim_body ~vars tg]: expects the target [tg] to point at a marked sequence.
     Then it will change all the declaraed variables inside that sequence  based on [vars]
     Either the user can give a list of variables together with their new names, or he can give the postifx
     that's going to be assigned to all the declared vairables. *)
(*========================*)
(** [bind ~fresh_name ~args tg]: expects the target [tg] to point at a function call,
    Then it will just call bind args and bind_intro.
    Basically this tranasformation just binds a variable to the targeted function call
    and its arguments.*)
(*========================*)
(** [inline ~resname ~vars ~args ~keep_res ~delete ~debug tg]: expects the target [ŧg] to point at a function call
    Then it will try to inline that function call. If it's possible this transformation tries to
    perform as many simplifications as possible.


    -------------------------------------------------------------------------------------------------------------
    This transformation handles the following cases:

    Case 1: Function call belongs to a variable declaration:
      Ex:
      int f(int x){
        return x + 1;
      }
      int a = 10;
      int b = f(a);
    Case 2: Function call belongs to a write operation
      Ex:
      int f(int x){
        return x + 1;
      }
      int a = 10;
      a = f(a);
    Case 3: Function call is does not return any value(a call to a function of void type)
      Ex:
        void f(int& x){
          x = x + 1;
        }
        int& a = 10;
        f(a);
    Case 4: Function call belongs to a for loop component
      Ex:
        int f(int x){
          return x + 1;
        }
        int a = 10;
        for(int i = f(a); i < 20; i ++){
          ..
        }
    -------------------------------------------------------------------------------------------------------------
      STEPS:

      Step 1(Only for case 1):
        Mark the instruction that contains the function call as "__inline_instruction"

      Step 2:
        Bind [resname] variable to the function call, if [resname] was not provided by the user then "__TEMP_Optitrust" is going to be
        used as a temporary variable.

      Step 3:
        Mark the function call for easier targeting in case it hasn't been marked by previous transformation.

      Step 4:
        Create a special mark for the inline body of the function [body_mark = "__TEMP_BODY" ^ (string_of_int i)] for easy targeting that inilin
        function body.

      Step 5:
        Call [Function_basic.inline] with target being the marked function call.
        Note: This step detaches the binded declaration.

      Step 6:
        Function arguments are encoded as const variables and that's different from the encodings of declared variables.
        This introduces wrong struct and array accesses. To fix this [Accesses.intro] with target being the marked body genereated
        from [Step 5] is called.

      Step 7:
        Integrates the sequence from [Step 6] to its surrouding sequence.
        To avoid name clashes this transformation renames all the variables defined inside that sequence by using the rule defined
        by [vars] variable.
        Note: It's the users responsibility to introduce a good renaming strategy(See [Variable_core.rename module]).

      Step 8:
        Recall [Step 5] detaches the binded declaration. This step tries to attach that variable declaration with the value returned
        by the function.

      Step 9:
        TODO: Add all the steps when function inline was completely debugged

      TODO: when the system of intermediate steps for combi transformations is available,
      apply it to generate the intermediate steps for the inlining example.

   EXAMPLE:
    int g(int x, int y, int z, int w) {
  int p = x + x + y + z + w;
  return p + p;
}
int main1() { // initial step : target on g(..)
  int u = 1, v = 2, w = 3;
  int t = f(g(h(4), u, m(v, 2), (w + 1)));
}
int main2() { // Function_basic.bind_intro
  int u = 1, v = 2, w = 3;
  int r = [mymark:](g(h(4), u, m(v, 2), (w + 1)));
  int t = f(r);
}
int main3() { // Function.bind_args ~[cMark mymark]
  int u = 1, v = 2, w = 3;
  int a = h(4);
  int b = m(v, 2);
  int r = [mymark:](g(a, u, b, (w + 1)));
  int t = f(r);
}
int main4() { // Function_basic.inline ~[cMark mymark]
              // The mark gets moved to the surrounding declaration
  int u = 1, v = 2, w = 3;
  mymark: int r; // same as before, only you remove the initialization term
  mybody: {
    int p = ((((a + a) + u) + b) + (w + 1));
    r = (p + p);
  }
  int t = f(r);
}
int main5() { // Function.elim_body ~[cMark mark]
  int u = 1, v = 2, w = 3;
  int a = h(4);
  int b = m(v, 2);
  mymark: int r;
  int p = ((((a + a) + u) + b) + (w + 1));
  r = (p + p);
  int t = f(r);
}
int main6() { // Variable_basic.init_attach
  int u = 1, v = 2, w = 3;
  int a = h(4);
  int b = m(v, 2);
  int p = ((((a + a) + u) + b) + (w + 1));
  int r = (p + p);
}

Other example in the case of return:

int h(int x) {
  if (x > 0)
    return -x;
  return x;
}
int f1() {
  int a = 3;
  int r = [mymark:]h(a);
  int s = r;
}
int f2() { // result of Funciton_basic.inline_cal
    // generate goto, generate label, don't call init_attach
  int a = 3
  [mymark:]int r;
  if (a > 0) {
    r = -a;
    goto _exit;
  }
  r = a;
  _exit:;
  int s = r;
} *)
(*========================*)
(** [inline_def]: like [inline], but with [tg] targeting the function definition.
   All function calls are inlined, with [delete = true] as default. *)
(*========================*)
(** [beta ~indepth tg]: expects the target [tg] to be pointing at a function call or a function declaration whose
     parent trm is a function call. If its the first case then it will just call Function_basic.beta.
     If its the second case then this transformation will just redirect the target to the parent function call
     and then call Function_basic.beta.
     [indepth]: if true it will apply the beta reduction to all the descendants of [tg].
     [body_mark]: mark left in case this transformation is used as an intermediate step of an another transformation.


     Note: If [tg] points to a function call then similar to Function_basic.beta, this transformation can be considered as an
     alias of Function_basic.inline. If that's not the case then transformation will do something as the similar to the following:
     int a = (void f(int x) {return x})(3) --> int a = 3;.*)
(*========================*)
(** [beta ~indepth tg]: applies beta-reduction on candidate function calls that appear
    either "exactly at" or "anywhere in depth" in the target [tg], depending on the value of ~indepth. *)
(*========================*)
(** [use_infix_ops ~tg_ops]: expects the target [tg] to be pointing at an instruction that can be converted to
     an infix form, for example x = x + 1 can be converted to x += 1,
    [indepth]: if true then it will check all the descendants of [t] if there are any write operations to be transformed
    [allow_identity]: if true it stops the transformation from failing when it finds nodes that can't be transformed.*)
(*========================*)
(** [uninline ~fxt tg]: expects the target [tg] to be pointing at an instruction that is similar to the first instruction
    of the body of the function declared in [fct]. Let nb be the number of instruction on the body of [fct]. The transformation
    will put the targeted instruction together with the following (nb -1) instructions into a sequence marked with a mark.
    Now the stage is ready for applying the basic version of uninline. After calling that transformation and assuming that
    everything went fine we can now eliminate the introduced sequence. The arg [with_for_loop] should be set to true if the
    original function declaration contains a for loop.*)
(*========================*)
(** [insert ~reparse decl tg]: expects the relative target [t] to point before or after an instruction,
     then it will insert the function declaration [decl] on that location.
     To integrate the new declaration with the current AST [reparse] should be set to true. *)

(* //////////////////////////////
   //////  Loop_basic.ml  ///////
   ////////////////////////////// *)

      (*========================*)
(** [color nb_colors i_color tg]: expects the target [tg] to point at a simple for  loop,
   let's say [for (int i = start; i < stop; i += step) { body } ].
   [nb_colors] - an expression denoting the number of colors (e.g., ["2"]),
   [index] - denotes a fresh name to use as index for iterating over colors.

   In case [step = 1]:
   {@c[for (int index = 0; index < nb_color; index++) {
      for (int i = index; i < stop; i += nb_color) { body }]}.

   In the general case, it produces:
   {@c[for (int index = 0; index < nb_color; index++) {
      for (int i = index*step; i < stop; i += step*nb_color) { body }]}. *)
(*========================*)
(** [tile tile_size index tg]: expects the target [tg] to point at a simple loop,
   say [for (int i = start; i < stop; i += step) { body } ].
   divides - denotes a flag to know if tile_size divides the size of the array or not
   [tile_size] - denotes the width of the tile (e.g., ["2"])
   [index] - denotes a fresh name to use as index for iterating over tiles.
   [bound] - can be one of
      - TileBoundMin: generates a constraint of the form  [i < min(X, bx+B)]
      - TileBoundAnd: generates a constraint of the form [i <  X && i < bx+B]
      - TileDivides: generates a constraint of the form [i < X], which is only true if B divides X

   It produces:
   {@c[for (int index = 0; index < stop; index += tile_size) {
      for (int i = index; i < min(X, bx+B); i++) { body }]}. *)
(*========================*)
(** [collapse]: expects the target [tg] to point at a simple loop nest:
    [for i in 0..Ni { for j in 0..Nj { b(i, j) } }]
    And collapses the loop nest, producing:
    [for k in 0..(Ni*Nj) { b(k / Nj, k % Nj) } ]

    Correct if Ni >= 0 and Nj >= 0.

    This is the opposite of [tile].

    LATER:
    - could generate binders for i, j
    - could generate i := PROJ1(Ni, Nj, k), j := PROJ2(Ni, Nj, k)
    - then can simplify MINDEX2(Ni, Nj, PROJ1(Ni, Nj, k), PROJ2(Ni, Nj, k))
      |--> MINDEX2PROJ(Ni, Nj, k) |-- elim_mops --> k
    - can be generalized to collapsing N loops and deriving PROJNM and MINDEXNPROJ
    *)
(*========================*)
(** [hoist x_step tg]: expects [tg] to point at a variable declaration inside a
    simple loop. Let's say for {int i ...} {
        int x; [tg]
        ...
        x = ..
      }
    The targeted declaration should be detached, then the transformation it's going to introduce
    an array declaration right before the for loop that contains the targeted declaration.
    The declared array will have name [name], type the same as the one targeted by [tg] and the size
    of the array it's going to be equal to the [loop_bound -1]. All the variable occurrences are
    going to be replaced with array accesses at index the index of the for loop.

    [x_step] - denotes the array name that is going to hoist all the values of the targeted variable
    for each index of the for loop. *)
(*========================*)
(** [fission_on_as_pair]: split loop [t] into two loops

    [index]: index of the splitting point
    [t]: ast of the loop
    *)
(*========================*)
(** [fission_on]: split loop [t] into two loops

    [index]: index of the splitting point
    [t]: ast of the loop
    *)
(*========================*)
(** [fission tg]: expects the target [tg] to point somewhere inside the body of the simple loop
   It splits the loop in two loops, the spliting point is trm matched by the relative target.

   @correctness: Reads in new second loop need to never depend on writes on
   first loop after index i. Writes in new second loop need to never overwrite
   writes in first loop after index i. *)
(*========================*)
(** [fusion]: expects the target [tg] to point at a loop that is followed by another loop with the same range (start, stop, step).
  Merges the two loops into a single one, sequencing the loop bodies into a new loop body:

  for (int i = start; i < stop; i += step) {
    body1
  }
  for (int i = start; i < stop; i += step) {
    body2
  }

  -->

  for (int i = start; i < stop; i += step) {
    body1;
    body2
  }
 *)
(*========================*)
(** [grid_enumerate index_and_bounds tg]: expects the target [tg] to point at a loop iterating over
    a grid. The grid can be of any dimension.
    Loop  [tg] then is transformed into nested loops
    where the number of nested loops is equal to the number of dimensions.
      [index_and_bounds] - is a list of pairs, where each pair denotes the index and the bound
        of the loop iterating over a specific dimension.
    Ex: Assume A = X * Y * Z, and [index_and_bounds] = [("x","X");("y","y");("z","Z")] and the result is

      for (int a = 0; a < A; a++){        for (int x = 0; x < X; x++){
        .......                       =>    for (int y = 0; y < Y; y++){
      }                                       for (int z = 0; z < Z, z++){
                                                int a = ((x * Y) + y)*Z + z
                                                ...
                                              }
                                            }
                                          } *)
(*========================*)
(** [unroll ~braces ~my_mark tg]: expects the target to point at a simple loop of the shape
    for (int i = a; i < a + C; i++) or for (int i = 0; i < C; i++)
      then it will move the instructions out of the loop by replacing
      the index i occurrence with a + j in and j in the second case where
      j is an integer in range from 0 to C.

    Assumption: Both a and C should be declared as constant variables. *)
(*========================*)
(** [move_out_on trm_index t]: moves an invariant instruction just before loop [t],
    [trm_index] - index of that instruction on its surrouding sequence (just checks that it is 0),
    [t] - ast of the for loop.
  *)
(*========================*)
(** [move_out tg]: expects the target [tg] to point at the first instruction inside the loop
    that is not dependent on the index of the loop or any local variable.
    Then it will move it outside the loop.

    {@c[
    for (i) {
      tg;
      middle-instrs
    }
    will become
    tg; // or: if (range not empty) tg;
    for (i) {
      middle-instrs
    }
    ]}

    Correctness check:
    (1) tg is idempotent/not self-interfering:
      tg uses resources RO(A) and Uninit(B), it must not use any full permission,
      it must not have any produce apart from the RW coming from B
    (2) Duplicating tg after middle-instrs would be redundant:
      in middle-instrs usage, resources A and B can only be used in RO mode.
      This is equivalent to say that there is no usage of A or B in Uninit or RW mode in middle-instr.
    (3) If the loop range can be empty, adding the extra tg instruction must not change the behaviour.
      Three methods to handle that:
      - Add an if on tg outside the loop
      - Prove that the loop range is never empty
      - All resources in B are uninit in the loop contract
*)
(*========================*)
(** [move_out_alloc ~empty_range tg]: same as [move_out], but supports moving out an allocation
    instruction together with its corresponding deallocation (that must be at the end of the loop).

    TODO: generalize [move_out] and [move_out_alloc] to moving out the any first/last group of instrs (I1 and I2) at once, if (I2; I1) is a no-op.
  *)
(*========================*)
(** [unswitch tg]:  expects the target [tg] to point at an if statement with a constant condition
     (not dependent on loop index or local variables) inside a loop.  Then it will take that
      if statment outside the loop.

   @correctness: requires that the loop is parallelizable *)
(*========================*)
(** [to_unit_steps index tg]: expects target [tg] to point at a for loop
    [index] - denotes the new index for the transformed loop
        by default is an empty string. The reason for that is to check if the user
        gave the name of the new index of not. If not then [index] = unit_index
        where index is the index of the targeted loop.

    Assumption:
      The targeted loop should be of the form:
        for (int i = a; i < b; i+=B){ s += i },
        and it assumes that B divides (b-a). It then
        transforms the targeted loop into the following form:
          for (int index = 0; index < ...; index++) {
            int i = (a + (j * B));
            s += i;
           } *)
(*========================*)
(** [scale_range ~factor ?index tg]: expects target [tg] to point at a for loop
    [index]. [factor] denotes the factor by which indices are multiplied. *)
(*========================*)
(** [fold ~direction index start stop step tg]: expects the target [tg] to point at the first instruction in a sequence
    and it assumes that the sequence containing the target [tg] is composed of a list of instructions which
    can be expressed into a single for loop with [index] [direction] [start] [nb_instructions] and [step] as loop
    components. *)
(*========================*)
(** [split_range nb cut tg]: expects the target [tg] to point at a simple loop
    then based on the arguments nb or cut it will split the loop into two loops. *)
(*========================*)
(** [shift_on index kind]: shifts a loop index to start from zero or by a given amount. *)
(*========================*)
(** [shift index kind]: shifts a loop index range according to [kind], using a new [index] name.

  validity: expressions used in shifting must be referentially transparent, other expressions can be bound before using [Sequence.insert].
  *)
(*========================*)
(** [extend_range]: extends the range of a loop on [lower] and/or [upper] bounds.
   The body of the loop is guarded by ifs statements, doing nothing on the extension points.

   For this to be correct, the loop bounds must be extended, not shrinked.
  *)
(*========================*)
(** [rename_index new_index]: renames the loop index variable *)
(*========================*)
(** [slide]: like [tile] but with the addition of a [step] parameter that controls how many iterations stand between the start of two tiles. Depending on [step] and [size], some iterations may be discarded or duplicated.
*)
(*========================*)
(** [delete_void]: deletes a loop with empty body. *)

(* //////////////////////////////
   ////////  Stencil.ml  ////////
   ////////////////////////////// *)

      (*========================*)
(** [may_slide]: slides a stencil that writes to [written] with outer loop at path [p], so that tiles of [sizes] values are produced by inner loops, within outer loops progressing by [steps].

  Keeps outer loop index names and uses [written] as suffix for inner loop index names.
  Returns the list of created inner loop index names.
  *)

(* //////////////////////////////
   ////////  Arrays.ml  /////////
   ////////////////////////////// *)

      (*========================*)
(** [inline_constant] expects the target [decl] to point at a constant array literal declaration, and resolves all accesses targeted by [tg], that must be at constant indices.
For every variable in non-constant indices, this transformation will attempt unrolling the corresponding for loop.
  *)
(*========================*)
(** [elim_constant] expects the target [tg] to point at a constant array literal declaration, and resolves all its accesses, that must be at constant indices. Then, eliminates the array declaration.
  *)

(* //////////////////////////////
   ////////  Record.ml  /////////
   ////////////////////////////// *)

      (*========================*)
(** [set_explicit tg]: an extension to [Record_basic.set_explicit](see Record_basic.ml), contrary to the basic
    on this transformation supports automatic variable declaration detachment.
    vect v = {0,0}; becomes vect v; v.x = 0; v.y = 0; *)
(*========================*)
(** [set_implicit tg]: an extension to [Record_basic.set_implicit](see Record_basic.ml), contrary to the basic one
     this one expects that the target [tg] matches all the write operations that can be converted to a single
     write operation. *)
(*========================*)
(** [rename_field field ~into tg]: this is a specialization of [Record_basic.rename_fields]
      when one wants to rename only one field of a Record. [field] is the current field name
      [into] is the new name that is going to replace all the occurrences of field in the context of
      the targeted typedef Record. *)
(*========================*)
(** [align_field align pattern tg]: expects the target [tg] to be pointing at a typedef struct definition,
   then it will align all the fields that match [pattern] with [align] size. *)

(* //////////////////////////////
   //////  Record_core.ml  //////
   ////////////////////////////// *)

      (*========================*)
(** [set_explicit_on t]: transforms an assigment into a list of field assignments,
     [t] - ast of the assignment. *)
(*========================*)
(** [set_implicit t]: transform a sequence with a list of explicit field assignments into a single assignment,
      [t] - ast of the sequence containing the assignments. *)
(*========================*)
(** [contains_field_access f t]: checks if [t] contains an access on field [f] *)
(*========================*)
(** [inline_struct_accesses x t]: changes all the occurrences of the struct accesses to a field into a field,
      [x] - the name of the field for which the transformation is applied,
      [t] - ast node located in the same level as the stract declaration or deeper. *)
(*========================*)
(** [inline_struct_initialization struct_name field_list field_index t]: changes all struct in struct initializations,
      [struct_name] - the type of the struct that is being inlined,
      [field_list] - a list of fields from the original type of the struct,
      [field_index] - index of the field in the outer struct,
      [t] - ast node located in the same level as the main struct declaration or deeper. *)
(*========================*)
(** [reveal_field_at field_to_reveal index t]: reveals field [field_to_reveal] on its typedef struct definition,
     update all the initializations and accesses according to this change,
      [field_to_reveal] - field that is going to be revealed,
      [index] - index of the struct declaration inside the sequence it belongs to,
      [t] - trm corresponding to a typedef struct definition. *)
(*========================*)
(** [fields_order]: the order should be provided as argument to the transformation [reorder_fields]. *)
(*========================*)
(** [compute_bijection order fl]: based on the [order] given, computes the bijection
    of the indices after applying that order. *)
(*========================*)
(** [reorder_fields_at order index t]: reorders the fields of the struct [t] based on [order],
     [order] - order based on which the fields will be reordered,
     [t] - ast of the typedef Record. *)
(*========================*)
(** [inline_struct_accesses name field t]: transforms a specific struct access into a variable occurrence,
    [name] - name of the variable to replace the struct access,
    [field] - struct accesses on this field are going to be replaced with [name],
    [t] - ast node located in the same level as the variable declaration. *)
(*========================*)
(** [to_variables_at index t]: changes a variable declaration of type typedef struct into a list
      of variable declarations with types inherited from the fields of the underlying type,
      [index] - index of the declaration inside the sequence it belongs to,
      [t] - ast of the surrounding sequence of the variable declarations. *)
(*========================*)
(** [Rename]: a module used for renaming the struct fields. *)
(*========================*)
(** [rename]: instantiation of module [Rename]. *)
(*========================*)
(** [rename_struct_accesses struct_name renam t]: renames all struct accesses based on [rename],
      [struct_name] - the constructed type whose fields are going to be renamed,
      [rename] - a type used to rename the struct fields,
      [t] - any node in the same level as the struct declaration.*)
(*========================*)
(** [rename_fields_at index rename t]: renames struct fields in the typedef struct definitions,
      [index] - the index of the struct declaration in the sequence [t],
      [rename] - a type used to rename the fields,
      [t] - the ast of the sequence which contains the struct declaration. *)
(*========================*)
(** [update_fields_type_aux pattern ty t]: changes the current type for all the struct fields,
      that are matched with [pattern] and whose type can be changed by [typ_update].
      [pattern] - regular expression to match struct fields,
      [typ_update] - function that modifies only specific types,
      [t] - the ast of the typedef definition. *)
(*========================*)
(** [simpl_proj_on t]: transforms all expression of the form {1, 2, 3}.f into the trm it projects to,
      [t] - ast of the node whose descendants can contain struct initialization list projections. *)
(*========================*)
(** [Struct_modif]: a module for defining struct modifications. *)
(*========================*)
(** [modif_accesses struct_name arg t]: modify struct accesses,
    [old_and_new_fields] - used to replace some specific fields,
    [struct_name] - used for checking the type of the struct access,
    [arg] - see Struct_modif module,
    [t] - trm corresponding to the surrounding sequence of the targeted typedef. *)
(*========================*)
(** [struct_modif_at new_fields f_get f_set use_annot_of index t],
     [arg] - Struct_modif type,
     [index] - index of the typedef on its surrounding sequence,
     [t] - ast of the main sequence containing the typedef definition. *)
(*========================*)
(** [change_field_access_kind_on acc_kind f t]: changes the access_kind for field [f] to [acc_kind] of class or struct [t]. *)
(*========================*)
(** [method_to_const_on method_name t]: converts the [method_name] method to a a const one,
    if the targeted method is already const than this transformation does nothing.
    [method_name] - the name of the method that's going to be converted.*)

(* //////////////////////////////
   /////  All_transfos.ml  //////
   ////////////////////////////// *)

      
(* //////////////////////////////
   ///////  Loop_swap.ml  ///////
   ////////////////////////////// *)

      (*========================*)
(** [swap_basic tg]: expects the target [tg] to point at a loop that contains an
   immediately-nested loop. The transformation swaps the two loops. *)
(*========================*)
(** [swap tg]: expects the target [tg] to point at a loop that contains an
  immediately-nested loop. The transformation swaps the two loops.
  Also handles ghosts that may be around the nested loop, and __sequentially_reads parallelization.

  {v
  for i: <<< tg @outer_loop_m
    __sreads(H);
    GHOST_BEGIN(gp, g);
    for j: <<< @inner_loop_m
      body;
    GHOST_END(gp);

  -- Resources.loop_parallelize_reads -->

  GHOST_BEGIN(hp, ro_fork_group(..));
  pfor i: << @outer_loop_m
    GHOST_BEGIN(gp, g);
    for j: <<< @inner_loop_m
      body;
    GHOST_END(gp);
  GHOST_END(hp);

  -- Ghost_pair.elim_all_pairs_at -->

  GHOST_BEGIN(hp, ro_fork_group(..));
  pfor i: << @outer_loop_m
    g();
    for j: <<< @inner_loop_m
      body;
    g_rev();
  GHOST_END(hp);

  -- Loop_basic.fission_basic -->

  GHOST_BEGIN(hp, ro_fork_group(..));
  pfor i:
    g();
  pfor i:
    for j: <<< @inner_loop_m
      body;
  pfor i:
    g_rev();
  GHOST_END(hp);

  -- Ghost.embed_loop -->

  GHOST_BEGIN(hp, ro_fork_group(..));
  ghost ([&] {
    pfor i:
      g();
  });
  pfor i:
    for j: <<< @inner_loop_m
      body;
  ghost ([&] {
    pfor i:
      g_rev();
  });
  GHOST_END(hp);

  -- Ghost_pair.reintro_pairs_at -->

  GHOST_BEGIN(hp, ro_fork_group(..));
  GHOST_BEGIN(gp, [&] {
    pfor i:
      g();
  }, [&] {
    pfor i:
      g_rev();
  });
  pfor i:
    for j: <<< @inner_loop_m
      body;
  GHOST_END(gp);
  GHOST_END(hp);

  -- Loop.swap_basic -->

  GHOST_BEGIN(hp, ro_fork_group(..));
  GHOST_BEGIN(gp, [&] {
    pfor i:
      g();
  }, [&] {
    pfor i:
      g_rev();
  });
  for j: <<< @inner_loop_m
   pfor i:
      body;
  GHOST_END(gp);
  GHOST_END(hp);
  v}

   *)

(* //////////////////////////////
   ////////  Typedef.ml  ////////
   ////////////////////////////// *)

      
(* //////////////////////////////
   /////////  Label.ml  /////////
   ////////////////////////////// *)

      
(* //////////////////////////////
   ///////  Accesses.ml  ////////
   ////////////////////////////// *)

      
(* //////////////////////////////
   //////  Align_core.ml  ///////
   ////////////////////////////// *)

      (*========================*)
(** [def_on vec_align t]: adds the alignas attribute to the  declaration [t],
      [vec_align] - alignment size,
      [t] - ast of the declaration. *)

(* //////////////////////////////
   //////  Specialize.ml  ///////
   ////////////////////////////// *)

      (*========================*)
(** [variable]: adds a specialized execution path when [var == value]. *)
(*========================*)
(** [variable_multi]: repeats [variable] transfo to create multiple specialized paths.

    - [mark_then]: function to create a mark based on [var, value]
    *)

(* //////////////////////////////
   /////  Sequence_core.ml  /////
   ////////////////////////////// *)

      (*========================*)
(** [insert_at index code t]: inserts trm [code] at index [index] in sequence [t],
    [index] - a valid index where the instruction can be added,
    [code] - instruction to be added as an arbitrary trm,
    [t] - ast of the outer sequence where the insertion will be performed. *)
(*========================*)
(** [delete_at index nb_instr t]: deletes a number of instructions inside the sequence starting
      from [index] and ending at ([index] + [nb]),
      [index] - starting index,
      [nb] - number of instructions to delete,
      [t] - ast of the outer sequence where the deletion is performed. *)
(*========================*)
(** [intro_at index nb t]: regroups instructions with indices falling in the range \[index, index + nb) into a sub-sequence,
       [mark] - mark to insert on the new sub-sequence,
       [label] - a label to insert on the new sub-sequence,
       [index] - index where the grouping is performed,
       [nb] - number of instructions to consider,
       [t] - ast of the outer sequence where the insertion is performed.

  Correct if the variables bound in the sub-sequence are only used locally, i.e. there is no scope interference with the outer sequence continuation.
*)
(*========================*)
(** [elim_on t]: inlines an inner sequence into the outer one,
      [t] - ast of the sequence to be removed.

  This function does not perform the inlining, but simply tags the inner sequence as 'nobrace'.

  Correct if there is no name conflict between the variables bound in the  inner sequence and the ones used in the outer sequence continuation, i.e. there is no scope interference with the outer sequence continuation. This is checked by the code responsible for removing 'nobrace' sequences.
  *)
(*========================*)
(** [wrap_on visible mark t]: surround [t] with a sequence,
    [mark] - mark to be added on the introduced sequence,
    [visible] - a flag on the visibility of the introduced sequence,
    [t] - any trm. *)
(*========================*)
(** [unwrap_on t]: the opposite of [wrap]
     [t] - a term that corresponds to a sequence with a single item in t. *)
(*========================*)
(** [split_at index t]: splits [t] into two sequences,
      [index] - the location where the splitting is done,
      [is_fun_body] - flag used when splitting function bodies,
      [t] - trm that corresponds to the the targeted sequence. *)
(*========================*)
(** [partition_on blocks braces]: partitions sequence [t] into a list of sequences,
      [blocks] -  a list of integers, where each integer denotes the size of the partition blocks,
      [braces] - flag on the visibility of the generated sequences,
      [t]: trm corresponding to a sequence. *)
(*========================*)
(** [shuffle_on braces t]: transposes a a list of partitioned sequences,
      [braces] - denotes a flag on the visibility of the added sequences,
      [t] - the ast of the complex sequence of blocks. *)

(* //////////////////////////////
   //////  Ghost_pair.ml  ///////
   ////////////////////////////// *)

      (*========================*)
(** [elim_all_pairs_at gen_mark p]: remove all the pairs inside the sequence pointed by the given path, generating marks to merge them back later.
  Returns the list of eliminated pair in their sequence order.
  In the returned list each ghost pair is represented by a triplet containing:
  - the variable name of the pair,
  - the generated mark on the ghost replacing the ghost_begin,
  - the generated mark on the ghost replace the ghost_end. *)

(* //////////////////////////////
   /////////  Loop.ml  //////////
   ////////////////////////////// *)

      (*========================*)
(** [rename]: instantiation of Rename module *)
(*========================*)
(** [hoist_alloc_loop_list]: this transformation is similar to [Loop_basic.hoist], but also supports undetached
   variable declarations, hoisting through multiple loops, and inlining array indexing code.
    [tmp_names] - pattern used to generate temporary names
    [name] - name of the hoisted matrix
    [inline] - inlines the array indexing code
    [loops] - loops to hoist through (expecting a perfect nest of simple loops),
              where [0] represents a loop for which no dimension should be created,
              and [1] represents a loop for which a dimension should be created.
  *)
(*========================*)
(** [hoist ~name ~array_size ~inline tg]: this transformation is similar to [Loop_basic.hoist] (see loop_basic.ml) except that this
    transformation supports also undetached declarations as well as hoisting through multiple loops.
    [inline] - inlines the array indexing code
    [nest_of] - number of loops to hoist through (expecting a perfect nest of simple loops)

    for (int l = 0; l < 5; l++) {
      for (int m = 0; m < 2; m++) {
        int x = ...;
      }
    }
    --> first hoist
    for (int l = 0; l < 5; l++) {
      int* x_step = MALLOC1(2, sizeof(int));
      for (int m = 0; m < 2; m++) {
        int& x = x_step[MINDEX1(2, m)];
        x = ...;
      }
    }
    --> second hoist
    int* x_step_bis = MALLOC2(5, 2, sizeof(int));
    for (int l = 0; l < 5; l++) {
      int*& x_step = x_step_bis[MINDEX2(5, 2, l, 0)];
      for (int m = 0; m < 2; m++) {
        int& x = x_step[MINDEX1(2, m)];
        x = ...;
      }
    }
    --> final
    int* x_step_bis = MALLOC2(5, 2, sizeof(int));
    for (int l = 0; l < 5; l++) {
      for (int m = 0; m < 2; m++) {
        int& x = x_step_bis[MINDEX2(5, 2, l, m)];
        x = ...;
      }
    }
 *)
(*========================*)
(** [hoist_instr_loop_list]: this transformation hoists an instructions outside of multiple loops using a combination of
  [Loop_basic.move_out], [Instr_basic.move], and [Loop.fission].
  [loops] - loops to hoist through (expecting a perfect nest of simple loops),
            where [0] represents a loop for which no dimension should be created,
            and [1] represents a loop for which a dimension should be created.
*)
(*========================*)
(** [hoist_decl_loop_list]: this transformation hoists a variable declaration outside of multiple loops
   using a combination of [hoist_alloc_loop_list] for the allocation and [hoist_instr_loop_list] for the initialization. *)
(*========================*)
(** [hoist_expr_loop_list]: this transformation hoists an expression outside of multiple loops
   using a combination of [Variable.bind] to create a variable and [hoist_decl_loop_list] to hoist the variable declaration. *)
(*========================*)
(** [hoist_expr]: same as [hoist_alloc_loop_list], but allows specifying
   loop indices that the expression does not depend on in [indep],
   and specifying where to hoist using [dest] target. *)
(*========================*)
(** [hoist_expr]: same as [hoist_expr_loop_list], but allows specifying
   loop indices that the expression does not depend on in [indep],
   and specifying where to hoist using [dest] target. *)
(*========================*)
(** [shift ~index kind ~inline]: shifts a loop index according to [kind].
- [inline] if true, inline the index shift in the loop body *)
(*========================*)
(** [extend_range]: like [Loop_basic.extend_range], plus arithmetic and conditional simplifications.
   *)
(*========================*)
(** [fusion nb tg]: expects the target [tg] to point at a for loop followed by one or more for loops.
    Merge them into a single loop.

    [nb] - denotes the number of sequenced loops to consider.
    [nest_of] - denotes the number of nested loops to consider.
    [adapt_fused_indices] - attempts to adapt the indices of fused loops using [Loop.extend_range] and [Loop.shift], otherwise by default the loops need to have the same range.
  *)
(*========================*)
(** [fusion_targets tg]: similar to [fusion] except that this transformation assumes that [tg] points to multiple
    not neccessarily consecutive for loops.
    All targeted loops must be in the same sequence.

  [into] - Specifies into which loop to fuse all other.
    Otherwise, fuse into the first loop in the sequence.
  [nest_of] - denotes the number of nested loops to consider.

  LATER ?(into_occ : int = 1)
  *)
(*========================*)
(** [move_out ~upto tg]: expects the target [tg] to point at an instruction inside a for loop,
    then it will move that instruction outside the for loop that it belongs to.
    In case of nested loops the user can specify the index of the upmost loop before which
    the instructions is going to be moved to.*)
(*========================*)
(** [move before after loop_to_move]: move one loop before or after another loop in
     a "sequence"(not in the context of Optitrust) of nested loops.
     [before] - a default argument given as empty string, if the user wants to move
     [loop_to_move]: before another loop then it should use this default argument with the
                     value the quoted loop index
     [after] - similar to [before] but now is the index of the loop after whom we want to move [loop_to_move]. *)
(*========================*)
(** [unroll tg] unrolls a loop nest and perform arithmetic simplification on the resulting trm. *)
(*========================*)
(** [reorder ~order tg]:  expects the target [tg] to point at the first loop included in the [order]
    list, then it will find all the nested loops starting from the targeted loop [tg] and
    reorder them based on [oder].

    Assumption:
      All loops have as bodies blocks of code(sequences).

    @correctness: correct if loops are parallelizable. *)
(*========================*)
(** [bring_down_loop]: given an instruction at path [p_instr], find a surrounding
   loop over [index] and bring it down to immediately surround the instruction.
   In order to swap imperfect loop nests, local variables will be hoisted ([Loop.hoist]),
   and surrounding instructions will be fissioned ([Loop.fission]).

   Returns a mark on the instruction at path [p_instr].
   *)
(*========================*)
(** [reorder_at ~order tg]: expects the target [tg] to point at an instruction that is surrounded
   by [length order] loops, and attempts to reorder these loops according to [order].
   The loops do not have to be perfectly nested. In order to swap imperfect loop nests,
   local variables will be hoisted ([Loop.hoist]),
   and surrounding instructions will be fissioned ([Loop.fission]).

   Example: order = [k; i; j]

  for i:
    ia;
    ib;
    for j:
      ja;
      for k:
        ka;
        kb; <--- tg @m
      jb;
    ic;

  --- bring_down_loop j @m --->

  for i:
    ia;
    ib;
    ...
    for j: <--- (loop_path 2)
      ja;
    ...
    for k: <--- @m_instr2
      for j:
        ka;
        kb; <--- @m / @m_instr
    ...
    for j:
      jb;
    ...
    ic;

  --- bring_down_loop i @m2 --->

  for i: <--- (loop_path 2)
    ia;
    ib;
    ...
    for j:
      ja;
  ...
  for k: <--- @m_instr2
    for i:
      for j: <--- @m2 / @m_instr
        ka;
        kb; <--- @m
  ...
  for i:
    for j:
      jb;
    ...
    ic;

  --- bring_down_loop k @m3 --->

  for i: <--- (loop_path 2)
    ia;
    ib;
    ...
    for j:
      ja;
  ...
  for k: <--- (loop_path)
    for i: <--- @m3 / @m_instr
      for j: <--- @m2
        ka;
        kb; <--- @m
  ...
  for i:
    for j:
      jb;
    ...
    ic;

  --- done ---

   *)
(*========================*)
(** [fold ~index ~start ~sstep ~nb_instr tg]: similar to [Loop_basic.fold] (see loop_basic.ml) except that
    this one doesn't ask the user to prepare the sequence of instructions. But asks for the first instructions and
    the number of consecutive instructions [nb_instr] that can be converted into a single loop.
   @correctness: always correct, as we can map all intermediate predicates
   to numbered predicates on the loop. *)
(*========================*)
(** [fold_instrs ~index ~start ~step tg]: similar to [fold] except that this one asks the user to provide a generic target
     that can match all the instructions that can be converted into a single loop. *)
(*========================*)
(** [unroll_first_iterations nb tg]: expects the target [tg] to be pointing at a simple loop;
   it extracts the sequences associated with the [nb] first iterations before loop.
   . *)
(*========================*)
(** [unroll_first_iteration tg]: expects the target [tg] to be pointing at a simple loop, it
   extracts the sequence associated with the first iteration before the loop. *)
(*========================*)
(** [unfold_bound tg]: inlines the bound of the targeted loop if that loop is a simple for loop and if that bound
    is a variable and not a complex expression. *)
(*========================*)
(** [grid_enumerate  ~indices tg]: similar to [Loop_basic.grid_enumerate](see loop_basic.ml) but this one computes
     the bounds automatically under the assumption that the bound of the targeted loop is given as a product of
    the bounds for each dimension. *)
(*========================*)
(** [change_iter iterator_function main_loop_function tg]:  TODO ARTHUR spec *)
(*========================*)
(** [collapse]: expects the target [tg] to point at a simple loop nest:
    [for i in 0..Ni { for j in 0..Nj { b(i, j) } }]
    And collapses the loop nest, producing:
    [for k in 0..(Ni*Nj) { b(k / Nj, k % Nj) } ]

    This is the opposite of [tile].
    *)
(*========================*)
(** [slide]: like [tile] but with the addition of a [step] parameter that controls how many iterations stand between the start of two tiles. Depending on [step] and [size], some iterations may be discarded or duplicated.
*)
(*========================*)
(** [slides]: like [slide], but operating on a nest of multiple loops and putting all loops over elements inside the bunch of loops over tiles. *)
(*========================*)
(** [delete_void]: deletes a loop nest with empty body.

  [nest_of] - number of perfectly nested loops to delete
  *)
(*========================*)
(** [delete_void]: deletes all loop nests with empty body. *)

(* //////////////////////////////
   /////////  Align.ml  /////////
   ////////////////////////////// *)

      (*========================*)
(** [alloc vec_align tg]: expects the target [tg] to point at a call at a OptiTrust MALLOC macro,
    then them will convert it to an aligned one with alignment size [vec_align]. *)

(* //////////////////////////////
   //////  Arith_core.ml  ///////
   ////////////////////////////// *)

      (*========================*)
(** [has_mark_nosimplf t]: check if [t] should be skipped by the simplifier or not.*)
(*========================*)
(** [transform aop inv pre_cast post_cast u t]: shifts or scale the right hand
    side of a set operation with term [u]
    [aop] - a flag to decide if the arithmetic operation should be Arith_scale
       or Arith_shift
    [inv] - a flag for the sign(plus or minus) of shifting
    [u] - shift size
    [pre_cast] - casting of type [pre_cast] performed on the right hand side of the
      set operation before shifting
    [post_cast] - casting of type [post_cast] performed after shifting
    [t] - the ast of the set operation *)
(*========================*)
(** [apply op arg t]: applies binary_operation [op] on [t] with the second  argument of the operation being [arg],
    [op] -the binary operation to apply.
    [arg] - the second operand of [op].
    [t] - the first argument in the performed operation. *)
(*========================*)
(** [expr]: expression type, it may be a literal expression, an atom expression or an arithmetic expression *)
(*========================*)
(** [exprs]: list of expressions *)
(*========================*)
(** [wexprs]: weighted list of expressions *)
(*========================*)
(** [wexpr]: weighted expression *)
(*========================*)
(** [Atom_map]: a map from atom ids to the corresponding terms *)
(*========================*)
(** [atom_map]: atom map for storing atoms *)
(*========================*)
(** [no_atoms]: empty atom map *)
(*========================*)
(** [expr_int n] produces the integer value [n] *)
(*========================*)
(** [expr_float f] produces the float value [n] *)
(*========================*)
(** [expr_one typ] produces either [expr_int 1] or [expr_float 1.0] depending on the type *)
(*========================*)
(** [expr_atom id] produces a variable [id], denoting an arbitrary subterm form ast.ml *)
(*========================*)
(** [expr_sum [(w1,e1);(w2,e2)]] produces [w1*e1 + w2*e2] *)
(*========================*)
(** [expr_sum_nonweighted es] produces [e1 + e2 + ... + en] *)
(*========================*)
(** [expr_neg e1] produces [-e1] *)
(*========================*)
(** [expr_add e1 e2] produces [e1 + e2] *)
(*========================*)
(** [expr_sub e1 e2] produces [e1 - e2] *)
(*========================*)
(** [expr_prod [(w1,e1);(w2,e2)]] produces [e1^w1 * e2^w2] *)
(*========================*)
(** [expr_prod_nonweighted es] produces [e1 * e2 * ... * en] *)
(*========================*)
(** [expr_pow e1 w] produces [e1 ^ w] *)
(*========================*)
(** [expr_mul e1 e2] produces [e1 * e2] *)
(*========================*)
(** [expr_div e1 e2] produces [e1 / e2];
   if arguments are integers, then [e1] is assumed to be divisible by [e2] *)
(*========================*)
(** [expr_binop op e1 e2] produces the operation [op e1 e2] *)
(*========================*)
(** [expr_div_floor e1 e2] produces the integer division [e1 / e2], rounded below *)
(*========================*)
(** [normalize_one e]:
   - collapses nested sums onto a single sum, and likewise for nested products
   - turns a product of an expression with a constant integer as a weighted
     expression in the parent sum
   - eliminates products and sums with a single expression of weight one
   - eliminates products and sums with an empty list
   - eliminates elements with weight zero
   - eliminates +0 is sums and *1 in produts
   - simplifies interger-division by 1
   - simplifies modulo operations applied to zero
   - simplifies binary shifting operations by zero
   *)
(*========================*)
(** [normalize e]: applies [normalize_one] in a bottom up fashion.
   Its definition appears further below. *)
(*========================*)
(** [is_one e]: checks if e == 1 *)
(*========================*)
(** [parens_if_neg n d]: if [n] is negative then it add parentheses around [d] *)
(*========================*)
(** [expr_to_string atoms e]: convert an expression to a string, in AST form *)
(*========================*)
(** [expr_to_math_string atoms e]: converts an expression to a string, using mathematical notations *)
(*========================*)
(** [identity] transformation *)
(*========================*)
(** [apply_bottom_up]: is a combinator that takes a transformation and applies it recursively,
   bottom up through a term. *)
(*========================*)
(** [normalize e]: applies [normalize_one] in a bottom up fashion *)
(*========================*)
(** [cleanup_true]: perform cleanup after transformation; else [cleanup_false] *)
(*========================*)
(** [recurse_true]: apply transformation in depth; else [recurse_false] *)
(*========================*)
(** [apply_bottom_up_if]: is a combinator for either applying a transformation recursively
   or applying it only at the top level, according to the [recurse] argument.
   If the [cleanup] argument is true, then after each call to the transformation,
   the operation [normalize_one] is called. *)
(*========================*)
(** [apply_bottom_up_debug e]: function used only for debugging purposes *)
(*========================*)
(** [apply_bottom_up_if_debug]: function used only for debugging purposes *)
(*========================*)
(** [create_or_reuse_atom_for_trm atoms t]: auxiliary function for [trm_to_naive_expr]*)
(*========================*)
(** [trm_to_naive_expr]: conversion of a trm from the AST into an expr, plus a map that for each atom gives
    the corresponding term *)
(*========================*)
(** [trm_to_expr t]: convert trm [t] to an expression*)
(*========================*)
(** [expr_to_trm atoms e]: converts expr [e] to trm  *)
(*========================*)
(** [cancel_div_floor_prod wes n] simplifies
   [Expr_div (Expr_prod wes) e] into [Expr_prod wes'].
   It returns [Some wes'], or [None] if no simplification is possible *)
(*========================*)
(** [gather_one e]: regroups similar expression that appear inside a same product
    or sum. For example, [2 * e1 + (-1)*e1] simplifies to [e1] and
    [e1 * e2 * e1^(-1)] simplifies to [e2].
    Also changes [(a / b) / c] into [a / (b * c)], and simplifies
    [(a1 * c * a2) / (b1 * c * b2)] into  [(a1 * a2) / (b1 * b2)]
    where [c] is a common item to the numerator and divisor (order-insensitive).
    This includes simplifications of [a / a] to [1] and of [(a*b)/a] to [b]. *)
(*========================*)
(** [gather_common recurse_bool e]: apply [gather_one] in a full expression
    if recurse is set to true *)
(*========================*)
(** [gather] and [gather_rec] can be passed as arguments to [Arith.simpl] *)
(*========================*)
(** [expand_one e]: expands a sum that appears inside a product.
    For example, [e1 * (e2 + e3)] becomes [e1 * e2 + e1 * e3].
    It can also expand e.g. [(e1 + e2) * (e3 + e4) * (e5 + e6)].
    The function performs nothing if no expansion can be performed.
    At the very end, it applies [normalize] to the result.*)
(*========================*)
(** [expand_common recurse e]: calls [expand_one] recursively, calling the [gather] operations
    after each step. *)
(*========================*)
(** [expand] and [expand_rec] can be passed as arguments to [Arith.simpl] *)
(*========================*)
(** [list_find_at_most_one f l] takes a list [l] and returns zero or one item
    satisfying [f], and returns the list of remaining items.
    LATER: implement as tail rec *)
(*========================*)
(** [euclidian] exploits (a/q)*q + (a%q) = a.
The transformation looks for sums that contain:
  (1) the item [a % q] with multiplicative factor 1
  (2) the item [a / q] with multiplicative factor [q].
      where the division is an integer division with flooring.
Then it replaces the two terms with [a].

NOTE: must check [a] nonnegative.

LATER: generalize to: [a%q] with multiplicative factor k

implementation: extract all the pairs (a,q) that appear as [a%q] with
multiplicity 1, then for each of them, either extract [a/q] with multiplicity [q]
and put back [a] in the sum, or put back the original [a%q] in the sum.
*)
(*========================*)
(** [wexpr_is_numeric (w,e)] returns true if [e] is a constant [Expr_int]
   or [Expr_float]. *)
(*========================*)
(** [wexpr_is_int (w,e)] returns true if [e] is a constant [Expr_int]. *)
(*========================*)
(** [compute_power_int n w] computes [n^w] *)
(*========================*)
(** [compute_power_double f w] computes [f^w] *)
(*========================*)
(** [compute_wexpr_sum wes] assumes all items in [wes] to satisfy [wexpr_is_numeric],
   and it returns a single item [(w,e)] describing the numerical result.
   It is either of the form [(n, Expr_int 1)] in case the result is the integer [n],
   or of the forme [(1, Expr_float f)] in case the result is the float value [f].
   When the sum involves only integers, the result is an integer;
   if, however, the sum involves at least one double, it is a double. *)
(*========================*)
(** [compute_wexpr_prod wes] is similar to [compute_wexpr_sum], but for products.
   It returns a single item [(w,e)] describing the numerical result.
   It is either of the form [(1, Expr_int n)] in case the result is the integer [n],
   or of the forme [(1, Expr_float f)] in case the result is the float value [f]. *)
(*========================*)
(** [compute_one e]: performs simplification of operations between known constants.
    For example, [4 + a + 3] becomes [a + 7].
    For example, [4.3 + a + 3] becomes [a + 7.3].
    For example, [10 / 2.5 - 3] becomes [1.0].
    The operation [normalize] is called at the end. *)
(*========================*)
(** [compute] can be passed as arguments to [Arith.simpl] *)
(*========================*)
(** [map_on_arith_nodes tr t]: applies arithmetic simplification [tr] in depth of [t]*)
(*========================*)
(** [simplify indepth f t]: converts node [t] to an expression, then applies the
     simplifier [f], then it converts it back to a trm
    params:
      [f]: simplifier function
      [t]: the node on which the simplifications should be performed
    return:
      update t with the simplified expressions
  LATER: should [simplify false f t] fail if [t] is not an application of prim_arith? *)
(*========================*)
(** [check_int_compare cmp t1 t2] tries to statically check that [cmp t1 t2] always holds. *)

(* //////////////////////////////
   ////  Variable_basic.ml  /////
   ////////////////////////////// *)

      (*========================*)
(** [fold ~at tg]: expects the target [tg] to point at a variable declaration,
      [at] - denotes a target where the folding is done. If empty the folding operation
             is performed on all the ast nodes in the same level as the
             declaration or deeper, by default [at] = []. *)
(*========================*)
(** [unfold ~mark ~at tg]: expects the target [tg] to be pointing at a variable declaration,
    then it will replace all the occurence of this variable inside the terms pointed by the
    target [~at] by the variable definition.
      [mark] - a mark to be added everywhere the variable was unfolded,
      [at] - denotes a target where the unfolding is done. If empty the operation
            is performed on all the ast nodes in the same level as the
            targeted declaration or deeper, by default [at] = [],
*)
(*========================*)
(** [inline ~mark ~at tg]: expects the target [tg] to be pointing at a variable declaration,
    then it will find all the occurrences of that variable and replace them with its initial value.
      [mark] - a mark to be added everywhere the variable was inlined,
      [delete_decl] - if true, then the declaration is removed during the inlining.
*)
(*========================*)
(** [rename ~into tg]: expects the target [tg] to be pointing at a declaration, then it will
    rename its declaration and all its occurrences. *)
(*========================*)
(** [init_detach tg]: expects the target [tg] to point at a variable initialization.
   It then splits the instruction into a variable declaration and a set operation. *)
(*========================*)
(** [init_attach const tg]: expects the target [tg] to point at a variable declaration,
    Then it will search inside the sequence which contains the variable declaration
    for an unique assigment. Then it will replace that assignment with a new initialized
    variable declaration.
    [const] -denotes a booleean to decide if the new declaration is constant or not. *)
(*========================*)
(** [local_name_on mark curr_var var_typ local_var t] declares a local
  variable [local_var] and replaces [curr_var] with [local_var] in [t].
    - [curr_var]: the replaced variable
    - [var_typ]: the type of the variable
    - [local_var]: the name of the local variable to be declared
    - [t]: the modified term that contains [curr_var].
  *)
(*========================*)
(** [local_name ~var var_typ ~local_var tg] declares a local
  variable [local_var] and replaces [var] with [local_var] in
  the term at target [tg].

  TODO: compose as Instr.insert alloc; Storage/Cell.reuse;
  *)
(*========================*)
(** [delocalize array_size neutral_element fold_operation tg]: expects the target [tg] to point to
    a block of code of the following form
      T a

   { T x = a; // mendatory format for first instruction

      for (int i = ...)
         x++;

      a = x;  // mendatory format for last instruction
   }@nobrace then
   Then it will transform it into:
       T a

   {
      { T x[N];
         x[0] = a;
         for (k = 1; k < N; k++)
            x[k] = 0;
      }@nobrace

      for (int i = ...)
         a++;

      { a = 0;
         for (k = 1; k < N; k++)
            a = a + x[k];  // could be a += if exists
         }@nobrace

   }@nobrace.
   [index] - denotes the index of the two added loops
   [array_size] - denotes the size of the array inside the block
   [ops] - the delocalize operation, it can be an arithmetic delocalization or an object delocalization
    of the array declared inside the block. *)
(*========================*)
(** [change_type new_type tg]: expects [tg] to point a variable declaration, then it will change the type of
    that variable with [new_type]. *)
(*========================*)
(** [insert ~constr ~name ~typ ~value tg]: expects the target [tg] to point at any relative location in a sequence
     then it will insert a variable declaration on that location,
     [const] - if true, then the inserted variable is going to be immutable, otherwise mutable,
     [reparse] - if true it will reparse the full ast after applying the trasnformation,
     [value] - initial value of the inserted variable,
     [name] - name of the inserted variable,
     [typ] - typ of the inserted variable;.

    NOTE: if initialization [value] is not provided then the declaration will be un-initialized. *)
(*========================*)
(** [subst ~subst ~space tg]]: expects the target [tg] to point at any trm that could contain an occurrence of the
    variable [subst], then it will check for occurrences of the variable [subst] and replace is with [put]. *)
(*========================*)
(** [elim_reuse]: given a targeted variable declaration [let x = get(y)], eliminates the variable
  declaration, reusing variable [y] instead of [x].
  This is correct if [y] is not used in the scope of [x], and can be uninit after the scope of [x].

  TODO: think about the relationship between this, [reuse], [elim_redundant], and [local_name].
  local_name should Instr.insert alloc; Storage.reuse; Storage.read_last_write; Instr.delete x = x
  *)
(*========================*)
(** [bind ~const ~mark fresh_name tg]: expects the target [tg] to be pointing at any trm, then it will insert a variable declaration
      with name [fresh_name] just before the instruction that contains the target [tg], and replace the targeted trm with an occurrence
      of the variable [fresh_name].
      [const] - if true the binded variable will be immutable, otherwise mutable,
      [mark_let] - an optional mark for the created instruction
      [mark_occ] - an optional mark for the introduced occurrences
      [mark_body] - mark used for marking the body of the targeted trm,
      [typ] - type of the binded variable, needed when the type can't be deducted from the targeted trm,
      [fresh_name] - name of the binded variable. *)
(*========================*)
(** [to_const tg]: expects the target [tg] to be point at a variable declaration, then it will search inside
      the same scope if there are any write operations on that variable.
      If that's the case then the tranformation will fail(for safety reasons).
      Otherwise, first switch the mutability of that variable and then replace all get operations on that variable with its intialization
      value.
  @correctness: always correct if the result type checks. *)
(*========================*)
(** [to_nonconst tg]: expects the target [tg] to be point at a variable declaration,
      If the variable is mutable then does nothing, otherwise change the mutability of the targeted variable to a mutable one,
      and replace all the variable occurrences with a get operation containing that occurrence. *)
(*========================*)
(** [simpl_deref ~indepth tg]: expects the target [tg] to be pointing at expressions of the form  *(&b), &( *b) in depth
    if [indepth] is set to true or at the give target if [indepth] is false.*)
(*========================*)
(** [exchange var1 var2 tg]: expects the target [tg] to point at an instruction that contains both the
    variable [var1] and [var2], then it will try to swap all the occurrences of [var1] with [var2]. *)
(*========================*)
(** [ref_to_pointer tg]: expects thee target [tg] to be pointing at a reference declaration, then it will convert
    this reference into a pointer. *)
(*========================*)
(** [ref_to_var tg]: expects the target [tg] to point at a refernce declaration,
     then it will convert it into a simple variable declaration. *)

(* //////////////////////////////
   /////////  Apac.ml  //////////
   ////////////////////////////// *)

      (*========================*)
(** [constify tg]: expects target [tg] to point at a function definition. It
   constifies function arguments and functions whenever it is possible depending
   on data accesses, aliases and dependencies. *)
(*========================*)
(** [unfold_function_calls tg]: expects target [tg] to point at a function
   definition. It moves all function calls under target [tg] out of variable
   declarations and nested function calls.

    Example:

          int a = f(g(2));

    becomes:

          int __var_1;
          __var_1 = g(2);
          int __var_2;
          __var_2 = f(__var_1);
          int a = __var_2;

    However:

          int a;
          a = f(g(2));

    becomes:

          int a;
          int __var_1;
          __var_1 = g(2);
          a = f(__var_1);

    as the call to 'f' is already dissociated from the declaration of 'a'. See
    also comments within the function.
*)
(*========================*)
(** [parallel_task_group tg]: expects target [tg] to point at a function
    definition.

    The first step of the transformation consists in replacing return statements
    by gotos. At the beginning of the process, the function's body is wrapped
    into a sequence to which a mark is assigned.
    See [Apac_basic.use_goto_for_return] for more details.

    In the second step, we put the marked sequence into an OpenMP task group.
    See [Apac_basic.task_group] for more details. *)
(*========================*)
(** [decl_cptrs]: hashtable storing available varaibles and whether they are
    arrays. *)
(*========================*)
(** [get_delete_task ptrs]: returns an OpenMP task that deletes the variable
    contained in [ptrs]. Each variable has an InOut dependency. *)
(*========================*)
(** [heapify_nested_seq tg]: expects target [tg] to point at a sequence. Then, it
    will:
    - send on the heap all variables declared in this scope,
    - dereference all use of them (add  '*'),
    - add task to delete these variables before 'return', 'break' and 'continue'
      statements when needed. *)
(*========================*)
(** [get_all_vars acc t]: returns the list of variables used in
    application [t]. *)
(*========================*)
(** [get_dep var ty]: returns the dep of the typ [ty]. *)
(*========================*)
(** [sync_with_taskwait tg]: expects target [tg] to point at a function
    definition. Then, it will add 'taskwait' OpenMP pragmas before loop and
    conditional branches. It also adds it at the end of loops. *)
(*========================*)
(** [fun_loc]: function's Unified Symbol Resolution *)
(*========================*)
(** [dep_info]: stores data of dependencies of tasks. *)
(*========================*)
(** [fun_args_deps]: hashtable storing dep_info for each function argument. *)
(*========================*)
(** [dep_infos]: list of dep_info and the corresponding variables. *)
(*========================*)
(** [vars_depth]: hashtable storing the pointer depth of variable and its name.
    The name of the variable is in the key and the value in order to use
    Apac_basic.get_vars_data_from_cptr_arith. *)
(*========================*)
(** [is_base_type ty]: checks whether [ty] is a base type or not. *)
(*========================*)
(** [is_dep_in ty]: checks whether the OpenMP dependency type of [ty] is In. It
    does not handle auto! *)
(*========================*)
(** [get_functions_args_deps tg]: expects target [tg] to point at the root. Then,
    it will return the dep_info of each function argument. *)
(*========================*)
(** [count_unop_get t]: returns the number of nested unary operations of [t]. *)
(*========================*)
(** [is_unary_mutation t]: checks if [t] is a primitive unary operaiton that
    mutates a variable. *)
(*========================*)
(** [get_unary_mutation_qvar t]: returns fully qualified name of a variable
    behind unary mutatation operator. *)
(*========================*)
(** [get_apps_deps vd fad t]: gets the list of dependencies of trm [t]. Expects
    the transformation [unfold_funcall] to be applied before. The returned list
    may have duplicates or different dependencies for the same variable, so it
    must be post-processed. *)
(*========================*)
(** [sort_dep_infos dis]: returns the the list of dependencies for In and InOut
    dependencies from dep_info list [dis]. *)
(*========================*)
(** [get_cpagma_from_dep_infos dis_in dis_inout]: returns the cpragma
    corresponding to a task.

    [dis_in] - list of In dependencies,
    [dis_inout] - list of InOut dependencies. *)
(*========================*)
(** [insert_tasks_naive fad]: expects target [tg] to point at a function
    definition. Then, it will insert task in the function's body.

    This is a naive implementation that adds a task on every
    instruction/application. *)

(* //////////////////////////////
   ///////  Cast_core.ml  ///////
   ////////////////////////////// *)

      (*========================*)
(** [insert_on ty t]: adds a type in front of [t] to cast its current type to [ty],
      [ty] - the type on which [t] is going to be casted to,
      [t] - any ast node that allows casting. *)

(* //////////////////////////////
   /////////  Instr.ml  /////////
   ////////////////////////////// *)

      (*========================*)
(** [read_last_write ~write tg]: similar to [Instr_basic.read_last_write] except that this transformation
     tries to find the last write instead of asking explicitly for a target to that write,
     [write_mark]: a mark used by the transformation [inline_last_write] to mark the instruction that needs to
                   be deleted
     [write]: optional explicit target to the last write instruction
     Note: This transformation will fail in case it doesn't find a write to to the address that it reads from. *)
(*========================*)
(** [inline_last_write ~write ~write_mark tg]: similar to [read_last_write] except that this one
    deletes the last write operation. *)
(*========================*)
(** [accumulate tg]: similar to [Instr_basic.accumulate], however this transformation requires the user to use a
    different target. For this transformation the user needs to provide a target to the first instruction
    and the number [nb] of consecutive instructions that can be accumulated into a single one. *)
(*========================*)
(** [accumulate_targets tg]: similar to [Instr_basic.accumulate], the main difference is that this transformation
     expects the target [tg] to point to multiple instructions instead of a sequence that consists of multiple
     instructions that can be reduced to a single one. Here the regrouping is done automatically. *)
(*========================*)
(** [move_in_seq ~dest tg] perform the same actions as {!Instr_basic.move},
   but move the instructions with the ghost pairs they need around them. *)
(*========================*)
(** [gather_targets ~dest tg]: expects the target [tg] to point to one or more instructions, than
    it will move this instructions just before the instruction targeted by [dest].

    DEPRECATED with this interface:
    It should be replaced by a more robust equivalent using relative immediate children targets.

    Note: This transformation assumes that [tg]  is going to be pointing at multiple instructions,
    hence no need to provide the occurrecen target [nbMulti] before the main target.

    TODO: add flags to choose:
    - failing if the instructions cannot be put strictly side-to-side
    - allowing ghosts to be moved automatically if there is a dependence
    - allowing non-ghosts to be moved
   *)
(*========================*)
(** [move dest tg]: move the instructions pointed by [tg] to destination [dest]
   Note: If checking validity, dest must be in the same sequence as the moved instructions.
*)
(*========================*)
(** [move_out tg]: moves the instruction targeted by [tg], just before its surrounding sequence. *)
(*========================*)
(** [move_out_of_fun tg]: moves the instruction targeted by [tg] just befor the toplevel declaration function
    that it belongs to. *)
(*========================*)
(** [set_atomic tg]: just an alias to Omp.atomic tg, please refer to omp_basic.ml  line 9 *)
(*========================*)
(** [unset_atomic ty]: the opposite of [set_atomic]. *)

(* //////////////////////////////
   ///////  Sequence.ml  ////////
   ////////////////////////////// *)

      (*========================*)
(** [intro ~start ~stop ~nb ~on ~mark ~visible]: this is a high level function for inserting a subsequnece
    inside another sequence.
     [start] - denotes the target for the starting point of the sub-sequence, it should be used in conjunction with
              [nb] or [end].
     [end] - denotes the target for the end point of the sub-sequnece, it should be used in conjunction with
             [nb] or [start].
    [nb] - in the case when the user does not give the target of  the end point of the sequence he can give as
      argument the number of instruction to include comming after [start] or [end]. If used with [start] the sign
      of [nb] should be poistive otherwise it should be negative.
    [on] - denotes a single target to be isolated inside the sub-sequence. When [on] is used all the other
      except mark and visible shoold be left empty. *)
(*========================*)
(** [intro_targets tg]: expects the target [tg] to point at one or more consecutive instuctions
      then it will introduce a sequence that contains those instructions. *)
(*========================*)
(** [apply ~start ~stop ~nb f]: invokes [f mark] where the [mark] is attached to a temporary sequence created
   by [Sequence.intro ~start ~stop ~nb]. This sequence is eliminated immediately afterwards. *)

(* //////////////////////////////
   //////  Label_basic.ml  //////
   ////////////////////////////// *)

      (*========================*)
(** [add label tg]: adds a C-label named [label] to the front of the terms
   matching the target [tg].
   Does nothing if [label = no_label].

   @correctness: always correct. *)
(*========================*)
(** [remove label tg]: removes a C-label named [label] matched by th target [tg]. *)

(* //////////////////////////////
   ///////  Variable.ml  ////////
   ////////////////////////////// *)

      (*========================*)
(** [Rename]: this type is used for variable renaming, the user can choose between renaming all the variables
    on one block, by giving the suffix to add or he can also give the list of variables to
    be renamed where the list should be a a list of string pairs ex. (current_name, new_name). *)
(*========================*)
(** [map_f]: applies function [f] inside the rename type *)
(*========================*)
(** [fold ~at ~nonconst tg]: similar to [Variable_basic.fold] but this one can be forced
     to non-const variables by setting [nonconst] flag to true.
    @correctness: The folded expression should have no observable side-effect.
    Moreover, the expression should produce the same value as when it was
    evaluated the first time.
    Exists r such that for initialization and all replacement points we have
    { H } expr { fun r' => [r' = r] * H } with H beeing the local invariant

    NOTE: if applied on non const variable, the variable must additionaly not have
    been mutated between the replacement point and its initialization. *)
(*========================*)
(** [insert_and_fold]: expects the target [tg] to point at a relative location, then it inserts a new variable
     declaration at that location. The new declared variable is [name] with [typ] and [value].
     This variable will be folded on all the ast nodes that come after the declared variable. *)
(*========================*)
(** [local_name ~var ~local_var tg] replaces occurences of [var] with a new variable [local_var] around the term targeted by [tg].

  TODO:
  - [~use_heap = false] flag ?
  - deal with _Uninit cases in pre/post where reading/writing value around term is not needed.
  *)
(*========================*)
(** [delocalize var ~into ~mark ~arr_size ~neutral_element fold_operation tg]:
    expects the target [tg] to point at a for loop. Then it will surround this loop with a @nobrace
    sequence. After that it will apply another transformation called local other name. Which as the name
    suggests it will declare a new variable inside the targeted block and replace the current one with t he new one.
    Finally a last instruction is added to save all the changes to the old variable. Now the stage is
    ready to apply delocalize which basically declares a new array oof size [array_size]. Then it will
    transform the old loop into a parallel loop. And finally a reduction is done to save the result into
    the old variable.

    Ex:
    int ANY(int maxValue){ return 0;}

    const int N = 2;
    typedef int T;

    int main(){
      T a;
      for (int i = 0; i < N; i++){
         a++;
      }
      int y = 0;
      return 0;
    }
    STEP1:
      Introduce a new local name x for variable a and replace all its occurrences inside
      the section of interest which can be any ast node. In this particular case it is
      a for loop.

    int main() {
      T a;
      /*@section_of_interest*/ T x = a;
      for (int i = 0; (i < N); i++) {
        x++;
      }
      a = x; /*section_of_interest@*/
      int y = 0;
      return 0;
    }
    STEP 2:
      Apply_ the basic delocalize transformation
     int main() {
        T a;
        /*@section_of_interest*/ T x[N];
        x[0] = a;
        for (int dl_k = 1; (dl_k < N); dl_k++) {
          x[dl_k] = 0;
        }
        for (int i = 0; (i < N); i++) {
          x[ANY(N)]++;
        }
        a = x[0];
        for (int dl_k = 1; (dl_k < N); dl_k++) {
          a += x[dl_k];
        } /*section_of_interest@*/
        int y = 0;
        return 0;
    } *)
(*========================*)
(** [delocalize ~var ~into ~index ~mark ~ops ~array_size ~intos tg]: it's a continuation to the [delocalize] transformation
    that will unroll all the introduced loops from the basic delocalize transformation and convert the newly declared array
    to a list of variables namely for each index on variable, this variables should be given by the user through the labelled
    argument [vars]. *)
(*========================*)
(** [intro_pattern_array ~pattern_aux_vars ~const ~pattern_vars ~pattern tg]: expects the target [tg] to be
     pointing to expressions of the form [pattern], then it will create an array of coefficients for each
    [pattern_vars] and replace the current coefficients with array accesses. *)
(*========================*)
(** [detach_if_needed tg]: expects the target [tg] to be pointing at a variable declaration, then it will
    check if that declaration was already initialized or not, if that's the case than it will deatch that
    declaration, otherwise no change is applied. *)
(*========================*)
(** [reuse ~reparse space tg] expects the target [tg] to be poiting to a variable declaration, then it will
    remove that declaration and replace all its occurrences with [space]

   @correctness: correct if the previous variable space was never read after the reuse point. *)
(*========================*)
(** [renames rename tg]: expects [tg] to point at a sequence.
    [rename] can be either ByList l where l denotes a list of pairs on which
    each pair consists the current variable and the one that is going to replace it.
    Or AddSuffix s, in this case all the variables declared inside the targeted sequence
     are going to be renamed by adding the suffix [s] at the end of its current name. *)
(*========================*)
(** [unfold ~accept_functions ~simple_deref ~delete tg] expects the target [tg] to be pointing at a variable declaration.
    If the variable has a struct type then a mark is created and passed as an argument to Variable_basic.unfold
      on the next step, otherwise no mark needs to be created.
      We consider the following cases:
      1) If the targeted variable is mutable then we try to make it immutable by calling Variable.to_const.
         WARNING: This step will fail in the case when there are any write operations on the targeted varibles.
      2) If the transformation didn't fail in the first step we are sure that we are trying to inline a const variable
         and we can call safely Variable_basic.unfold
      3) If the targeted variable is a struct type then call Record_basic.simpl_proj to remove all the occurrences
          of struct initialization list, by projecting them on the field they are accessed.
          Ex: int v = {0,1} if we had v.x then Variable_basic.inline will transform it to {0, 1}.x which is non valid C code.
          After calling Record_basic.simpl_proj {0, 1}.x becomes 0 .
          Finally, if simple_deref is set to true then we will seach for all the occurrences of *& and &* and simplify them. *)
(*========================*)
(** [inline ~accept_functions ~simpl_deref tg]: similar to [unfold] except that this transformation
     deletes the targeted declaration by default. *)
(*========================*)
(** [inline_and_rename]: expects the target [tg] to point at a variable declaration with an initial value
    being another variable. Then it will inline the targeted variable. And rename variable that was the
    initialization value to the one we inlined

    Assumption:
      if the target [tg] points to the following instruction int y = x; then
      no occurrence of x appears after that instruction *)
(*========================*)
(** [elim_redundant ~source tg]: expects the target [tg] to be point at a variable declaration with an initial value being
    the same as the variable declaration where [source] points to. Then it will fold the variable at [source] into
    the variable declaration [tg] and inline the declaration in [tg]. *)
(*========================*)
(** [insert ~constr name typ value tg]: expects the target [tg] to point at a location in a sequence then it wil insert a
    new variable declaration with name: [name] and type:[typ] and initialization value [value].
    This transformation is basically the same as the basic one except that this has a default value for the type argument. *)
(*========================*)
(** [insert_list ~const names typ values tg]: expects the target [tg] to be poiting to a location in a sequence
    then it wil insert a new variable declaration with [name], [typ] and initialization [value] *)
(*========================*)
(** [insert_list_same_type typ name_vals tg]: inserts a list of variables with type [typ] and name and value give as argument in [name_vals]. *)
(*========================*)
(** [bind_multi ?all ?dest name tg] performs common subexpression elimination.
   Minimalistic example: [ f(e); g(e) ] becomes [ int x = e; f(x); g(x) ].
   It takes a target that aims at one or more expressions.
   If multiple expressions are targeted, the occurrences must correspond to equal expressions.
   LATER: If [~all:true] is provided, then we automatically consider replacing all occurrences
   within the syntactic scope of the new binding (and not just the ones targeted by [tg]).
   LATER: It takes as argument a target [dest] where to put the binding; if [dest] is not provided,
   then the point just before the instruction containing the first target is used.
   LATER: document arguments passed to [bind].
   LATER: add arguments for marks that could remain at the end
*)
(*========================*)
(** [bind_syntactic]: performs common subexpression elimination.
   Minimalistic example: [ f(a, b); g(a, b) ] becomes [ int x = a; int y = b; f(x, y); g(x, y) ].
   It takes a target that aims at one or more expressions.
   If multiple expressions are targeted, targeted expressions that are syntactically equal must correspond to semantically equal expressions.
   *)

(* //////////////////////////////
   //////  Arrays_core.ml  //////
   ////////////////////////////// *)

      (*========================*)
(** [inline_array_access array_var new_vars t]: change all the occurences of the array to variables,
    params:
      [array_var]: array_variable  to apply changes on
      [new_vars]: a list of variables, the variables at index i replaces and occurence of [array_var[i]]
      [t]: ast node located in the same level or deeper as the array declaration
    return:
        updated ast with the replaced array accesses to variable references. *)
(*========================*)
(** [to_variables_at new_vars t]: transform an array declaration into a list of variable declarations
      the list of variables should be entered by the user. The number of variables should correspond to
      the size of the arrys. The variable at index i in [new_vars] will replace the array occurrence
      at index i.
      (new_vars] - a list of strings of length equal to the size of the array,
      [index] - index of the instruction inside the sequence,
      [t] - ast of the surrounding sequence of the array declaration. *)
(*========================*)
(** [apply_tiling base_type block_name b x]: changes all the occurences of the array to the tiled form,
      [base_type] - type of the array
      [block_name] - new name for the array
      [b] - the size of the tile
      [x] - typvar
      [t] - ast node located in the same level or deeper as the array declaration

    assumptions:
    - if x is ty*, each array of type x is allocated through a custom function:
      x a = my_alloc(nb_elements, size_element)
    - x is not used in function definitions, but only in var declarations
    - for now: in any case, the number of elements is divisible by b. *)
(*========================*)
(** [tile_at block_name block_size index t]: transform an array declaration from a normal shape into a tiled one,
    then call apply_tiling to change all the array occurrences into the correct form.
      [block_name] - the name of the arrays representing one tile,
      [block_size] - the size of the tile,
      [index] - the index of the instruction inside the sequence,
      [t] - ast of the outer sequence containing the array declaration.

   Ex:
    t[N] -> t[N/B][B]  where t has the targeted type
    t[i] -> t[i/B][i%B] *)
(*========================*)
(** [apply_swapping x t]: swaps all array accesses that have type [t],
      [x] - typvar,
      [t] - an ast node which on the same level as the array declaration or deeper. *)
(*========================*)
(** [swap_at index t]: swap the dimensions of an array declaration,
     [index] - index of the array declaration on its surrouding sequence,
     [t] - AST of the surrouding sequence of the targeted array declaration. *)
(*========================*)
(** [aos_to_soa_rec struct_name sz t] : transforms an array of structures to a structure of arrays *)
(*========================*)
(** [detach_init_on t]: transform an initialized array declaration into a list of write operations
      for each one of its cells
    [t] - the array declaration. *)

(* //////////////////////////////
   /////  Record_basic.ml  //////
   ////////////////////////////// *)

      (*========================*)
(** [set_explicit tg]: expects the target [tg] to point at a set instruction where one struct
    instance has been assigned another struct instance. *)
(*========================*)
(** [set_implicit tg]: expects the target [tg] to point at a sequence containing
      a list of struct set assignments. And transforms it into a single struct assignment.
      So it is the inverse of set_explicit. *)
(*========================*)
(** [reorder_fields order tg]: expects the target to be pointing at typedef struct or class.
      then it changes the order of the fields based on [order].
      [order] - can be one of the following
        Move_before (x,fl) -> all the fields that belong to [fl] are moved before the field [x].
        Move_after (x,fl) -> all the fields that belong to [fl] are moved after the field [x].
        Reorder_all [fl] -> all the fields are reorder an will appear like [fl].

    @correctness: Correct if pointer arithmetic to field is replaced everywhere,
      might be impossible to prove in case of casts between types.*)
(*========================*)
(** [reveal_field ~reparse field_to_reveal_field tg]: expects the target [tg] to point at a typedef struct,
    then it will find [field_to_reveal_field] and it's underlying type and it will
    replace [field_to_reveal_field] with a list of fields rename comming from its underlying type. *)
(*========================*)
(** [reveal_fields fields_to_reveal_field tg]: an extension to the reveal_field transformation, this one
     is applied on multiple struct fields. *)
(*========================*)
(** [to_variables tg]: expects the target [tg] to point at a variable declaration of type typedef Record.
    Then it will transform this declaration into a list of variable declarations where the type
    of these variables is inherited from the type of the struct definition. All the struct_accesses
    are going to be changed to variable occurrences. *)
(*========================*)
(** [rename_fields rename tg] expects the target [tg] to point at a struct declaration,
    then it will rename all the fields that are matched when applying the type [rename]
    which can be a function to rename all the struct fields or only those that
    are matched by the [pattern] given as argument when the function [only_for] is used (see struc_core.ml). *)
(*========================*)
(** [applyto_fields_type ~reparse pattern typ_update tg]: expects the target [tg] to point at a
    struct definition, then it will update all the struct field types whose identifier matches [pattern]. *)
(*========================*)
(** [update_fields_type pattern ty tg]: expects the target [tg] to point at a struct declaration,
    then it will change the current type to [ty] for all the fields that are matched by [pattern]. *)
(*========================*)
(** [simpl_proj tg]: expects the target [tg] to point at any node whose descendants can contain struct
    initialization list projections. *)
(*========================*)
(** [struct_modif new_fields f_get f_set use_annot_of tg]: expects the target [tg] to point at a typedef struct,
    then it will replace its current fields with [new_fields]. After modifying the fields it will search for
    accesses of the targeted struct and modify them, if they are surrounded by a set operation it will apply
    [f_set] on that access otherwise [f_get] is going to be applied. *)
(*========================*)
(** [change_field_access_kind acc_kind f tg]: expects the target [tg] to point a typedef, then it will find
    field [f] at change its access kind to [acc_kind]. *)
(*========================*)
(** [make_all_members_public tg]: expects the target [tg] to point at a typedef struct or class.
    then it will transform all its members to public. *)
(*========================*)
(** [method_to_const method_name]: expects the target [ŧg] to be pointing at a typedef record definition.
    Then it will check if the method of that record definition is already a const method or not.
    If it's a const method then this transformation does nothing, otherwise it will transform that method to a const one.
    Note: If [method_name] is not specified by the user all the methods will be converted to const methods.*)

(* //////////////////////////////
   ///////  Omp_core.ml  ////////
   ////////////////////////////// *)

      
(* //////////////////////////////
   /////  Function_core.ml  /////
   ////////////////////////////// *)

      (*========================*)
(** [bind_intro_at my_mark index fresh_name vonst p_local t]: bind the variable [fresh_name] to the targeted function call,
      [my_mark] - put a mark on the targeted function call,
      [index] - index of the instruction that contains the targeted function call on its surrouding sequence,
      [const] - flag for the mutability of the binded variable,
      [p_local] - path from the instruction containing the function call to the function call itself,
      [t] - ast of the sequence that contains the targeted function call. *)
(*========================*)
(** [inline_at index body_mark p_local t]: inline a function call,
      [index] - index of the instruction containing the function call,
      [body_mark] - mark usef for the transflated body of the function,
      [p_local] - path from the instructions that contains the function call to the function call itself,
      [t] - ast of the sequence containing the function call. *)
(*========================*)
(** [use_infix_ops_on allow_identity t]: transforms an explicit write operation to an implicit one
      [allow_identity] - if true then the transformation will never fail
      [t] - ast of the write operation *)
(*========================*)
(** [uninline_on fct_decl t]: takes a function declaration [fct_decl], for example
   [void gtwice(int x) { g(x, x); }], and expects a term [t] that matches the body
   of the function, for example [g(3,3)]. It performs some matching to resolve [x]
   and returns the term [gtwice(3)], which is equivalent to [t] up to inlining. *)
(*========================*)
(** [trm_var_assoc_list to_map al]: creates a map from an association list wher keys are variables and values are trms *)
(*========================*)
(** [rename_args_on vl t]: renames arguments of function [t] and replace all the occurrences of its
    arguments of the args inside its body with the new names provided as arguments,
      [vl] - new arguments, can be [dummy_var] to avoid renaming.
      [t] - ast of the function declaration whose arguments are going to be altered. *)
(*========================*)
(** [replace_with_change_args_on new_fun_name arg_mapper t]: change the name of the called function and its arguments
      [new_fun_name] - the new name that is going to replace the current one,
      [arg_mapper] - a function to change the arguments. *)
(*========================*)
(** [dsp_def_at index arg func t]: changes the destination pasing style,
     [index] - index of the targeted function definition on its surrounding sequence,
     [arg] - the new argument to be added on the new function definition,
     [func] - name of the newly added function definition,
     [t] - ast of the original function definition. *)
(*========================*)
(** [dsp_call_on dps t]: changes a write operation with lhs a function call to a function call,
    [dsp] - the name of the function call, possibly empty to use the default name
    [t] - ast of the write operation. *)
(*========================*)
(** [get_prototype t]: returns the return type of the function and the types of all its arguments.*)

(* //////////////////////////////
   /////////  Cast.ml  //////////
   ////////////////////////////// *)

      
(* //////////////////////////////
   ////  Specialize_core.ml  ////
   ////////////////////////////// *)

      (*========================*)
(** [any_on e t]: replaces the function call [t] with [e]
      [e] - the expression replacing the call to function [ANY],
      [t] - ast of a call to function [ANY]. *)
(*========================*)
(** [choose_on selelct_arg t]: replaces the function call [t] with one of its arguments that satisfies
     the predicate  [select_arg],
      [select_arg] - a predicate on the index of the argument that should be choosed,
      [t] - ast of the call to function choose. *)
(*========================*)
(** [fun_def_aux spec_name spec_args t]: inserts a copy of the function definition [t], specializing
      one of its arguments based on the list [spec_args].
      [spec_name] - the name of the copy
      [spec_args] - an optional list of trms, telling the transformation which argss it shoudl specialize,
      [t] - ast of the function definition. *)
(*========================*)
(** [fun_call_on spec_name args_to_choose t]: replaces the function call [t] with another function call.
    [spec_name] - the name of the function that appears on the new call,
    [args_to_choose] - a list of booleans telling the transformation which of the current arguments
        from the current call should be kept.*)

(* //////////////////////////////
   /////////  Ghost.ml  /////////
   ////////////////////////////// *)

      
(* //////////////////////////////
   //////  Arith_basic.ml  //////
   ////////////////////////////// *)

      (*========================*)
(** [shift ~reparse ~inv ~pre_cast ~post_cast u tg]:  expects the target [tg]
    to point at a trm on which an arithmetic operation can be applied, then
    depending on the value of [inv] it will add or substract [u] to that trm.*)
(*========================*)
(** [scale ~inv ~pre_cast ~post_cast u] *)
(*========================*)
(** [apply op arg] expects the target [tg] to be pointing at any node of the ast
      then it applies the binary operation [op] at that node with the second argument
      of that operation being [arg] *)
(*========================*)
(** [simpl f] applies a arithmetic rewriting method from the module Arith_core:
   - gather  for grouping and cancelling out similar expressions in sums and produts
   - expand  for expanding products involving sums. *)
(*========================*)
(** [simpl_rec f tg] just an alias for simpl ~indepth:true tg *)
(*========================*)
(** [compose fs] returns the function obtained as the composition
    of the functions [fs] *)
(*========================*)
(** [simpls fs tg] is like simpl with the composition of the functions fs *)
(*========================*)
(** [simpls_rec f tg] just an alias for simpl ~indepth:true tg *)
(*========================*)
(** [simplify ~indepth tg] applies simpl with the operation being gathering of
    arithmetic experssions *)
(*========================*)
(** [clear_nosimpl tg]: clears all the marks on all the instructions that where
    skipped by the simplifier *)
(*========================*)
(** [nosimplf tg]: mark all the instructions targeted by [tg] as "__arith_core_nosimpl" *)
(*========================*)
(** [with_nosimpl tg f]: after marking all the nodes targeted by [tg] with mark "__arith_core_with_nosimpl", applies the
    transformation [f] on all the nodes matched  by [tg], after the transformation has been applied succesfully,
    it will clean all the introduced marks *)

(* //////////////////////////////
   ////////  If_core.ml  ////////
   ////////////////////////////// *)

      (*========================*)
(** [may_merge]: if [t] corresponds to two nested ifs, merge them using '&&'.
   *)

(* //////////////////////////////
   ////  Function_basic.ml  /////
   ////////////////////////////// *)

      (*========================*)
(** [delete tg]: delete the targeted function definition.
   Correct if the function is never used.

   Currently checked by verifying that the targets correspond to
   function definitions, and by retychecking the code *)
(*========================*)
(** [bind_intro ~fresh_name ~const ~my_mark tg]: expects the target [t] to point at a function call.
     Then it will generate a new variable declaration named as [fresh_name] with type being the same
     as the one of the function call, and initialized to the function call itself.
     If [const] is true the the bound variable will be declared as an immutable variable otherwise immutable.
     Then it will fold the newly declared variable.

     @correctness: correct if the new order of evaluation of expressions is
      not changed or does not matter. *)
(*========================*)
(** [inline ~body_mark tg]: expects the target [tg] to point at a function call inside a declaration
    or inside a sequence in case the function is of void type. Example:
          int r = g(a);
      or  g(a);

    Then it will replace that instruction with a nobrace sequence which is a sequence visible on the ast level.
    This sequence will be marked with [body_mark] and it will contain the body of the declaration of the called
    function targeted with [tg].
    Transformation steps:
       1) generate in that sequence the binding "int r", in case it is needed
          (if the original instructions featured a "int r = ..")

       2) replacing the name of the arguments with the expressions that were
           provided to the call.

          - Instructions of the form "return foo;" should be translated into
              "r = foo; goto __exit_body;"
          - Instructions of the form "return;" should be translated into
              "goto __exit_body;"
          - The "goto" is not needed in case the instruction is the "last one"
              of the body. To keep track of this, the recursive traversal function
              maintains a boolean flag "islast". This flag is true initially, but
              becomes false as soon as one enters the branch of a Trm_seq that is
            not the last branch. (see examples further below).
            - You can use a reference to save the information of whether at least
              one "goto" operation was generated.

        3) generate the exit label ("__exit_" ^ label) in case we observed the need
          for a goto during the translation of the body

     Example:

      int g(int x) {
        int y = x + x;
        return y + y;
      }
      int r = [target:]g(a)

     this result is:

        @nobrace{
          int r;
          body: {
             int y = a + a;
             r = y + y;
          }
        }

   @correctness: always works, and also needs to instantiate variables in the
   local invariants in the body. *)
(*========================*)
(** [beta ~body_mark tg]: similar to [function_inline] the main difference is that [beta] is used in the cases
    when the decaration of the function call can be founded at the targeted function call contrary to [inline]
    which will need to find first the toplevel declaration.  *)
(*========================*)
(** [use_infix_ops_at tg]: expects the target [tg] to point at an explicit set operation of the form x = x (op) a,
    then it will transform that instruction into x (op)= a. Ex: x = x + 1 --> x += 1. *)
(*========================*)
(** [uninline ~fct tg] expects the target [ŧg] to be pointing at a labelled sequence similar to what Function_basic.inline generates
    Then it will replace that sequence with a call to the fuction with declaration targeted by [fct]. *)
(*========================*)
(** [rename_args new_args tg]: expects the target [tg] to point at a function declaration, then it will rename the args of
     that function. If there are local variables declared inside the body of the function that have the same name as one
     of the function args then it will skip those variables on all their occurrences. *)
(*========================*)
(** [replace_with_change_args new_fun_name arg_mapper tg]: expects the target [tg] to point at a function call, then it will
    replace the name of the called function with [new_fun_name] and apply [arrg_mapper] to its arguments. *)
(*========================*)
(** [dsp_def ~arg ~func tg]: expects the target [tg] to point at a function definition, then it will
     inserts a new version of that definition whose return type is void.
    [arg] - is the name of the argument that's going to be inserted,
    [func] - the name of the new function that's going to be inserted. *)
(*========================*)
(** [dsp_call ~dsp tg]: expects the target [tg] to point at a function call whose parent trm is a write operation
    then it will convert that write operation into a function call.
    Let's say that the targeted function call is r = f(x, y);
    If [dsp] is the empty string, then "f_dsp" will be used as a name based on the original name "f".
    Note: This transformation assumes that dsp_def has been already applied to the definition of the called function. *)

(* //////////////////////////////
   ////////  Rewrite.ml  ////////
   ////////////////////////////// *)

      
(* //////////////////////////////
   //////////  If.ml  ///////////
   ////////////////////////////// *)

      
(* //////////////////////////////
   ////////  Reduce.ml  /////////
   ////////////////////////////// *)

      (*========================*)
(** [reduce]: constructs a call to [reduce]. *)
(*========================*)
(** [reduce_inv]: returns the list of arguments of a call to [reduce]. *)
(*========================*)
(** [elim_basic tg]: eliminates a call to [reduce], expanding it to a for loop. *)
(*========================*)
(** [elim_inline tg]: eliminates a call to [reduce], expanding it to an inlined expression.
    TODO: later, implement this as combi (1. unroll; 2. inline accumulator; 3. simplify zero add)
    *)
(*========================*)
(** [elim_basic tg]: eliminates a call to [reduce], expanding it to a for loop.

  - [unroll]: whether the reduction loop should be unrolled
  - [inline]: whether the reduction variable should be inlined (implies [unroll])
  *)
(*========================*)
(** [slide_basic tg]: given a target to a call to [set(p, reduce)] within a perfectly nested loop:
    [for i in 0..n { set(p, reduce(... i ...)) }]
    allocates a variable outside the loop to compute next values based on previous values:
    [alloc s = reduce(... 0 ...); set(p[i := 0], s); for i in 1..n { set(s, f(s)); set(p, s) }]
  *)
(*========================*)
(** [slide tg]: given a target to a call to [set(p, reduce)] within a perfectly nested loop:
    [for i in 0..n { set(p, reduce(... i ...)) }]
    allocates a variable outside the loop to compute next values based on previous values:
    [alloc s = reduce(... 0 ...); for i in 1..n { set(s, f(s)); set(p, s) }]

    TODO: generate check that n > 0
  *)

(* //////////////////////////////
   //////////  Omp.ml  //////////
   ////////////////////////////// *)

      (*========================*)
(** [ensure_header ()]: insert omp.h header at top of the file *)
(*========================*)
(** [set_num_threads threadnum tg]: sometimes omp_get_num_threads doesn't work with gcc for that we need to apply the following trick
     #pragma omp parallel
     {
       #pragma omp single
       nbThreas = omp_get_num_threads();
     } *)
(*========================*)
(** [parallel_for ~clause ~collapse tg]: when collapse is provided as argument then
     clause will not be taken into account*)

(* //////////////////////////////
   /////  Rewrite_basic.ml  /////
   ////////////////////////////// *)

      (*========================*)
(** [equiv_at rule]: expects the target [tg] to point at a trm where [rule] can be applied to transform
   that trm into a similar one defined by the rule itself. *)
(*========================*)
(** [compute tg]: expects the target [tg] to point at an arithmetic operation then it will try to simlplify it. *)
(*========================*)
(** [compute_inside tg]: expects the target [tg] to point at any trm in the ast that could contain some inner
   arithmetic operations, then it will try to simplify them by calling [compute] on that trm. *)

(* //////////////////////////////
   ////  Accesses_basic.ml  /////
   ////////////////////////////// *)

      (*========================*)
(** [transform ~reparse f_get f_set tg]: expects the target [tg] to point at a trm inside a set or a get operation.
    Then the transformation will search for the first get or set operation surrounding the targeted trm and call
    the transform core transformation on that trm. If the first founded operation was a get operation then [f_get]
    will be applied on the node represented by target [tg]. If it was a set operation then [f_set] will be applied
    on the second argument of the targeted node. *)
(*========================*)
(** [scale ~inv ~factor tg]: this transformation just calls the [transform] function  with [f_get] and [f_set] args
   defined as a multiplication and a division operation respectively. If [inv] is set to true then these two
   operations will be swapped. *)
(*========================*)
(** [shift ~inv ~factor tg]: this transformation just calls the [transform] function with [f_get] and [f_set] args
   defined as a multiplication and a division respectively. If [inv] is set to true then these two operations
   will be swapped. *)
(*========================*)
(** [intro tg]: expects the target [tg] to be pointing at any node that could contain struct accesses, preferably
   a sequence, then it will transform all the encodings of the form struct_get (get (t), f) to
   get (struct_access (t, f)) . *)

(* //////////////////////////////
   ///////  If_basic.ml  ////////
   ////////////////////////////// *)

      (*========================*)
(** [insert_aux single_branch index cond t]: takes one or two instructions and create an if statement or an if else
     statment when [single_brnach] is true,
      [cond] - condition of the if statement given as string as a trm,
      [t] - ast of the outer sequence containing the instruction. *)
(*========================*)
(** [insert cond tg]: expects the target [tg] to point at an instruction
    inside a sequence. Then it will create an if statement with the condition entered by the user
    and both its "then" and "else" branches will contain the same instruction.
    [cond] - the code which will appear as the condition in the if statement.
    [no_else_branch] - if true then the inserted if statement will not contain an else branch.
    Note:
      If [cond] is given as arbitrary string the flag [reparse] should be set to true. *)

(* //////////////////////////////
   /////////  Marks.ml  /////////
   ////////////////////////////// *)

      (*========================*)
(** [add_fake_instr m tg]: adds a mark [m] at the location of the target [tg] as a fake sequence item.
   NOTE: if m = "" then does nothing. *)
(*========================*)
(** [remove_fake_instr m tg]: remove a fake sequence item mark. *)

(* //////////////////////////////
   /////  Typedef_basic.ml  /////
   ////////////////////////////// *)

      (*========================*)
(** [fold ~at tg]: expects the target [tg] to point at a typedef, then it will fold its definition,
    [at] - denotes a target where folding is done. If empty
           the folding operation is performed on all the ast nodes in the same
           level as the declaration or deeper. *)
(*========================*)
(** [unfold ~delete ~at tg]: expects the target [tg] to point at a typedef, then it inlines its definition,
    [delete] - denotes a flag for telling if the declaration should be kept or no,
    [at] - denotes a target where inlining is done,
           if empty the inlining operation is performed on all the ast nodes in the
           same level as the declaration or deeper, by default [at] = []. *)
(*========================*)
(** [insert_copy name tg]: expects the target [tg] to point at a typedef, then copies the content
      of the body of typedef at gives to it the name [name]. *)
(*========================*)
(** [insert name td_body]: expects target [tg] to point at a relative location inside a sequence
    then it will insert a typedef declaration on that location.
    [name] - is the new type name while
    [td_body] - is the kind of typedef we are going to declare.
                It can be an alias a product(for struct declarations), a sum type or an enum. *)

(* //////////////////////////////
   ///////  Expr_core.ml  ///////
   ////////////////////////////// *)

      (*========================*)
(** [replace_fun_aux name t]: changes the current function call to another function call where the name has
      been changed to [name],
      [name] - name of the function replacing the targeted one,
      [t] - ast of the function call trm. *)
(*========================*)
(** [view_subterms_on stringreprs ro]: prints the string representations of all the subterms of [t]?
   assumes terms to carry [Annot_stringreprid], using the ids from the table [stringreprs].
   See [Expr_basic.view_subterms] for details on how this is achieved. *)

(* //////////////////////////////
   ///////  Apac_core.ml  ///////
   ////////////////////////////// *)

      (*========================*)
(** [lvar] Labelled variable type. We make use of this type to make a difference
   between class member variables. Indeed, in the case of a class member, the
   associated variable is represented always by [this]. As a result every member
   of a given class ends up with [this] as name and the same identifier (which
   changes only from one class method to another).

   For example, let us consider:

   class C {
     int a; int b;
     int f() { int c = a / b; return c; }
     int g() { return a * b; }
   }

   In [f], both [a] and [b] will have [this] as name as well as the same
   identifier, 75 for example. In [g], the identifier will change, let us say to
   42, but then both [a] and [b] in [g] will have [this] as name and 42 as
   identifier. The only way to make a difference between two class memebers
   within a member method is to look at the labels (strings) associated with
   each occurence. To take into account this situation within the constification
   process, we use a variable type extended with a label.
*)
(*========================*)
(** [lvar_to_string] returns a string representation of the labelled variable
   [lv]. *)
(*========================*)
(** [LVar] and [LVar_Hashtbl]: specific type of hash tables where the keys are of
   type [lvar]. *)
(*========================*)
(** [const_arg]: an argument constification record. *)
(*========================*)
(** [const_fun]: a function constification record. *)
(*========================*)
(** [const_funs]: type for a hash table of [const_fun]. The keys are functions
   represented by terms of type [var]. *)
(*========================*)
(** [const_aliases]: type for hash table of argument aliases.

   Pointers and references used within a function definition might be linked to
   the same data as the arguments of the function, i.e. they represent aliases
   to the function's arguments. Therefore, when we constify an argument, we must
   constify its aliase(s) too. To keep trace of argument aliases, we use a hash
   table data type where the key is the [lvar] of the alias and the value is a
   pair of a [lvar] and an [int]. The [lvar] element corresponds to the function
   argument being aliased. The [int] element gives the pointer degree of the
   alias, if the latter is a pointer, e.g. the pointer degree of [int ** tab] is
   2. *)
(*========================*)
(** [const_unconst]: type of stack of function arguments, except for objects,
   that must not be constified.

   At the beginning of the constification process, we assume that every function
   argument in every function can be constified, see
   [build_constification_records]. Then, we perform an analysis of dependencies
   between function arguments and write operations involving the latter within
   the body of corresponding functions, see [identify_mutables]. When the
   analysis concludes that a function argument is written to, it shall modify
   the associated constification record so as to mark the argument as
   non-consitifiable. The same goes also for class methods modifying sibling
   class members. Such methods should not be consitified either. See [const_arg]
   and [const_fun] for more details.

   However, the constification records are not modified directly during the
   analysis. Instead, we use a stack keeping trace of functions arguments and
   functions that shall be unconstified once the analysis terminates, see
   [to_unconst] below. The elements of the stack are pairs of the [var]
   identifying the target function and the [var] identifying its argument to
   unconstify. If the latter represents the [lvar] of the function, i.e. when
   both of the elements of the pair are equal [lvars], it means that the
   function itself should be unconstified. See
   [unconstify_mutables.unconstify_to_unconst_objects]. *)
(*========================*)
(** [const_unconst_objects]: type of stack of function arguments, represented by
   objects, that must not be constified.

   When an object is passed as argument to a function, it must not be constified
   if it is used to call a non-const class method which may modify one or more
   members of the class. However, this information is not know before the
   elements in a [const_unconst] stack are processed (see
   [unconstify_mutables]). This is why a second phase is necessary to process
   the elements of a [const_unconst_objects] stack (see [unconstify_mutables]).
   Here, an element is represented by a triplet of [var]s where the first [var]
   identifies the target function, the second [var] identifies the argument of
   that function to unconstify and the third [var] represents the class member
   method called on that argument. Then, during the unconstification process
   (see [unconstify_mutables]), if the method represented by the third [var]
   gets unconstified, the argument (second [var]) of the function targeted by
   the first [var] will be unconstified too. *)
(*========================*)
(** [const_mult]: type of hash table for multiple variable declarations that must
   be transformed into sequences of simple variable declarations while
   constifying one or more of the variables being declared.

   When a function argument is constified, it is necessary to propagate the
   constification to the aliases of the argument as well. This is done with the
   [Apac_basic.constify_arguments] function. This process is straightforward in
   the case of simple variable declarations. However, if an argument alias is
   declared within a multiple variable declaration, we have two different
   situations to consider:

   1) the inner type of the multiple variable declaration is already constant:

      // [i] is a constant integer, [j] is a mutable pointer to a constant
      // integer
      int const i, * j;

   2) the inner type of the multiple variable declaration is not constant:

      // [i] is a mutable integer, [j] is a mutable pointer to a mutable integer
      int i, * j;

   If we want to constify [j] in 1), it is possible without breaking the parent
   multiple variable declaration:

   int const i, * const j;

   In 2) we cannot fully constify [j], i.e. [int const * const j], without
   constifying [i], which is not always desirable. In this case, we have to
   split the parent multiple variable declaration into a sequence of simple
   variable declarations and then constify the declaration of [j] only. However,
   this transformation cannot be directly applied within
   [Apac_basic.constify_aliases_on]. Indeed, the introduction of a new sequence
   of instructions changes the scope of the declared variables. This can be
   prevented, but not in a term to term transformation function such as
   [Apac_basic.constify_aliases_on]. It must be done in a transformation
   function modyfing the AST by side-effect, i.e. a target to unit
   transformation, in which we can call [Nobrace_transfo.remove_after] to
   effectively remove the braces from the sequence of simple variable
   declarations so as to preserve their initial scope. Therefore,
   [Apac_basic.constify_aliases_on] should only mark those multiple variable
   declarations which need to be split into simple variable declarations and
   determine which declarations should be constified. To keep track of this
   information, so we can actually perform the transformations (using
   [Apac_basic.unfold_let_mult]), we use a hash table where the keys are the
   marks (strings) added to the concerned multiple variable declarations and
   where the values are lists of booleans. The latter have as many elements as
   there are declarations in the target multiple variable declaration. In other
   terms, a boolean value is associated to each variable declaration. If the
   value is [true], it means that the corresponding variable declaration should
   be constified. *)
(*========================*)
(** [const_mult]: hash table for multiple variable declarations that must be
   transformed into sequences of simple variable declarations while constifying
   one or more of the variables being declared. See [const_mult]. *)
(*========================*)
(** [typ_is_alias ty]: checks if [ty] is a user-defined alias to a basic type. *)
(*========================*)
(** [typ_get_alias ty]: if [ty] is a user-defined alias to a basic type, it
   returns the latter. *)
(*========================*)
(** [typ_get_degree ty]: computes and returns the pointer degree of the type
   [ty]. For example, for [int ** a], it returns 2. *)
(*========================*)
(** [trm_strip_accesses_and_references_and_get_lvar t]: strips [*t, &t, ...]
   recursively and if [t] is a variable, it returns the associated labelled
   variable. *)
(*========================*)
(** [typ_constify ty]: constifies [ty] by applying the 'const' keyword wherever
   it is possible. *)
(*========================*)
(** [trm_resolve_binop_lval_and_get_with_deref] tries to resolve the variable
   behind an lvalue and check whether it has been dereferences, i.e. following
   an array access or the use of [*]. Upon success, it returns the corresponding
   labelled variable. See [LVar] for more details on labelled variables. *)
(*========================*)
(** [trm_resolve_var_in_unop_or_array_access_and_get t] tries to resolve the
   variable (including the associated label if we are dealing with a class
   member variable) involved in a unary operation (++, --, & or get) or array
   access [t] and return it in a form of a labelled variable. *)
(*========================*)
(** [trm_resolve_pointer_and_aliased_variable t aliases]: tries to resolve
   pointer operation [t] and checks in [aliases] whether the resulting pointer
   is an argument or an alias to an argument. If the pointer operation succeedes
   and if the resulting pointer is an argument or an alias to an argument, the
   function returns the labelled variable corresponding to the resulting pointer
   as well as the aliased labelled variable. *)
(*========================*)
(** [trm_can_resolve_pointer t]: tries to resolve operation [t] to unique
   variable and returns [true] on success and [false] otherwise. *)
(*========================*)
(** [trm_let_update_aliases ?reference tv ti aliases]: checks in [aliases]
   whether the variable declaration specified by the typed variable [tv] and the
   initializaion term [ti] creates an alias to an already existing variable or
   function argument. If it is the case, it updates the alias hash table
   [aliases] accordingly and returns [1] if [tv] is a reference and [2] if [tv]
   is a pointer. Otherwise, it does nothing and returns [0].

   Note that when the function is applied on the elements of a simple variable
   declaration, i.e. a [trm_let], the [is_reference] function used to check
   whether the new variable is a reference has no effect. This is due to
   differences in representing simple [trm_let] and multiple [trm_let_mult]
   variable declarations in OptiTrust. In this case, the optional [reference]
   parameter can be used to ensure the test evaluates correctly. See a usage
   example in [const_compte_one]. *)
(*========================*)
(** [find_parent_typedef_record p]: goes back up the path [p] and returns the
   first term corresponding to a class or a structure definition, if any. We use
   this function to determine the parent class of a structure or a function in
   order to access to the member variables of that class or structure. *)
(*========================*)
(** [unconstify_mutables] proceeds in two phases, i.e. [unconstify_to_unconst]
   and [unconstify_to_unconst_objects]. At first, it pops out the elements from
   [to_unconst] and then from [to_unconst_objects] (global variables, see Part
   II.1) one by one and propagates the unconstification through the concerned
   function constification records in [const_records] (global variable, see Part
   II.1). *)
(*========================*)
(** [build_constification_records_on]: see [build_constification_records]. *)
(*========================*)
(** [build_constification_records]: expects the target [tg] to point at a
   function definition. It adds a new entry into [const_records] (global
   variable, see Part II.1) based on the information about the function. *)
(*========================*)
(** [identify_mutables_on t]: see [identify_mutables]. *)
(*========================*)
(** [identify_mutables tg]: expects the target [tg] to point at a function
   definition. It recurses over the body of the target function definition in
   order to figure out which function arguments and possible aliases to
   arguments should not be constified and adds them to the [to_unconst] stack
   (global variable, see Part II.1). *)

(* //////////////////////////////
   /////  Matrix_basic.ml  //////
   ////////////////////////////// *)

      (*========================*)
(** [intro_calloc tg]: expects the target [tg] to point at  a call to funciton alloc then it will
    replace this call with a call to CALLOC. *)
(*========================*)
(** [intro_malloc tg]: expects the target [tg] to point at a call to the function MALLOC,
      then it will replace this call with a call to MALLOC. *)
(*========================*)
(** [intro_mindex dim tg]. expects the target [tg] to point at an array access
    then it will replace that access to let say index i with an access at
    MINDEX (dim,i). *)
(*========================*)
(** [reorder_dims order tg]: expects the target [tg] to point at a call to ALLOC or MINDEX functions,
      then it will reorder their args based on [order], where [order] is a list of indices which the
      current args should follow. *)
(*========================*)
(** [insert_alloc_dim new_dim]: expects the target [tg] to point at call to ALLOC functions, then it will
      add a new arg at the begining of the list of args in the targeted call. *)
(*========================*)
(** [insert_access_dim new_dim new_index tg]: expects the target [tg] to point at an array access, then it will
    add two new args([new_dim] and [new_index]) in the call to MINDEX function inside that array access. *)
(*========================*)
(** [biject fun_name tg]: expectes the target [tg] to point at a function call, then it replaces the name
     of the called function with [fun_name]. *)
(*========================*)
(** [local_name ~mark var into tg]: expects the target to point at an instruction that contains
      an occurrence of [var] then it will define a matrix [into] whose dimensions will be the same
      as the one of [var]. Then we copy the contents of the matrix [var] into [into] and finally we
      free up the memory. *)
(*========================*)
(** [local_name_tile ?mark_accesses ?indices ~alloc_instr ?ret_var ~local_var tile tg]
  expects [alloc_instr] to point to an allocation of [var] and [tg] to point to a term [t]
  that contains occurrences of [var], and instead of [t]:
  + defines a matrix [local_var] that will correspond to a [tile] of matrix [var];
  + copies the contents from [var] to [local_var];
  + performs [t] where all accesses to [var] are replaced with accesses to [local_var];
  + copies the contents from [local_var] to [var];
  + frees up the memory of [local_var].

  - [mark_accesses] allows marking the produced [local_var] accesses.
  - [indices] allows providing explicit index variable names for the created for loop.
  - [ret_var] will receive the value [var].
  - [uninit_pre]/[uninit_post] specify whether [var] is uninit before/after [t], eliminating copies.
    For [uninit_pre = false], resource checks will fail if [var] is uninit before [t].
    For [uninit_pre = true], resource checks will fail if [var] is not uninit in [t],
    because [local_var] will be uninit.
    For [uninit_post = false], resource checks will always succeed.
    For [uninit_post = true], we check that [trm_ghost_forget_init ..var..] can be added after [t].

  TODO?
  - what if [alloc_instr] is not available? make it optional, retrieve var, dims, elem_ty, ...
    from MINDEX calls, as done by {! Matrix_stack_copy}
  - stack_alloc / heap alloc
  - check consistent API with {! Variable.local_name}
  - factorize and update {! Matrix_basic.local_name} with no tile
  - factorize with {! Matrix_stack_copy}
  *)
(*========================*)
(** [delocalize_aux dim init_zero acc_in_place acc any_mark labels index]: TODO  *)
(*========================*)
(** [delocalize ~init_zero ~acc_in_place ~acc ~dim ~index ~ops] a generalized version of variable_delocalize. *)
(*========================*)
(** [simpl_index_add]: simplifies an MINDEX(..) + MINDEX(..) expression,
   into a single MINDEX(..) expression, if the dimensions are compatible:

   MINDEX{N}  (            n1, .., nN,                 i1, .., iN) +
   MINDEX{N+M}(m1, .., mM,     .., m{N+M}, j1, .., jM,     .., j{N+M})
    = [if n{i} = m{i+M}]
   MINDEX{N+M}(m1, .., mM,     .., m{N+M}, j1, .., jM, i1 + j{M+1}, .., iN + j{N+M})

   For correctness, size and index expressions must be pure.
   *)
(*========================*)
(** [simpl_access_of_access]: simplifies &((&p[i0])[i1]) into &p[i0 + i1]

   TODO: should this be in another file?
   *)
(*========================*)
(** [intro_malloc0]: given a target to a sequence with a declaration allocating
   variable [x] on the stack, changes the declaration to use a MALLOC0 heap
   allocation, and adds an instruction to free the memory after all uses of
   [x] in the sequence.

   {
    T* x = new T*;
    ... uses x ...
    ...
   }
   -->
   {
     T* x = MALLOC0(sizeof(T));
     ... uses &x[MINDEX0()] ...
     free(x);
     ...
   }

   LATER: deal with control flow

   *)
(*========================*)
(** [stack_copy ~var ~copy_var ~copy_dims tg] expects [tg] to points at a term [t]
    that contains occurences of matrix [var], and instead of [t]:
    + defines a matrix [copy_var] that will correspond to the right-most [copy_dims]
      of [var], which are contiguous in memory;
    + copies the contents from [var] to [copy_var] using a memcpy;
    + performs [t] where all accesses to [var] are replaced with accesses to [copy_var];
    + copies the contents from [copy_var] to [var] using a memcpy;

  TODO:
  - [uninit_pre]/[uninit_post] as in {! Matrix.local_name_tile}

  This should be equivalent to using {! Matrix.local_name_tile}, converting its heap
  allocation to a stack allocation with something like {! Matrix.to_array},
  and converting copy loops to batch memory copies.
  *)
(*========================*)
(** [elim_mindex] expects target [tg] to point at a call to MINDEX,
  and replaces it with the flattened index computation.

   Equivalent to:
   Rewrite.equiv_at ~ctx:true "int d1, d2, i1, i2; ==> MINDEX2(d1, d2, i1, i2) == (i1 * d2 + i2)" tg
   Rewrite.equiv_at ~ctx:true "int d1, d2, d3, i1, i2, i3; ==> MINDEX3(d1, d2, d3, i1, i2, i3) == (i1 * d2 * d3 + i2 * d3 + i3)" tg
   [...]
   *)
(*========================*)
(** [storage_folding] expects target [tg] to point at a sequence defining matrix
   [var], and folds the [dim]-th dimension so that every index [i] into this matrix dimension is mapped to index [i % n].

   assumes that [i >= 0].
   *)
(*========================*)
(** [delete] expects target [tg] to point to a sequence defining matrix [var], and deletes it.
  Both allocation and de-allocation instructions are deleted.
  [var] should not be used anywhere, this is checked through var ids.
  TODO: additionnal check/invariant: this assumes that alloc/free dimensions are read-only expressions
   *)
(*========================*)
(** [read_last_write]: expects the target [tg] to pint at a matrix read operation, and replaces it with the value that was last written to this matrix index. The [write] target must correspond to this last write.
  For correctness, if [V] was written at index [i], reading [V[j/i]] should be equivalent to reading at index [j].
   *)

(* //////////////////////////////
   //////  Apac_basic.ml  ///////
   ////////////////////////////// *)

      (*========================*)
(** [mark_taskification_candidates_on]: see [mark_taskification_candidates]. *)
(*========================*)
(** [mark_taskification_candidates]: expects the target [tg] to point at a
   function definition. It marks the functions to taskify, i.e. the functions to
   the body of which we shall insert tasks. For now, we use a naive algorithm to
   determine the candidates for taskification. We consider every function
   performing at least two function calls. *)
(*========================*)
(** [task_group_on ~master t]: see [task_group] *)
(*========================*)
(** [task_group ~master t]: puts the instruction sequence of a function's body
    into an OpenMP task group, i.e. a block of instructions delimited by curly
    brackets and prepended with the OpenMP pragma '#pragma omp taskgroup'.

    Example:

          int g() {
            int a;
            f();
            return 0;
          }

          becomes

          int g() {
          #pragma omp taskgroup
          {
            int a;
            f()
          }
            return 0;
          }

    If [master] is true, the transformation creates a task group that will be
    executed only by one thread, the master thread, using the following pragmas:

      #pragma omp parallel
      #pragma omp master
      #pragma omp taskgroup

    Example:

      int main() {
        int a;
        f();
        return 0;
      }

      becomes

      int main() {
      #pragma omp parallel
      #pragma omp master
      #pragma omp taskgroup
      {
        int a;
        f()
      }
        return 0;
      }

    [master] - decides whether extra pragmas should be added for limiting the
               execution of the task group to only thread, the master thread;
    [t] - AST of a function body.
*)
(*========================*)
(** [use_goto_for_return_on mark t]: see [use_goto_for_return]. *)
(*========================*)
(** [use_goto_for_return mark]: expects the target [tg] to point at a function
    definition. It replaces potentially multiple return statements by a single
    return statement at the end of the function definition through the usage of
    gotos.

    First of all, the transformation wraps the function's body into a sequence
    and marks it with [mark] if [mark] <> "". Then,

    if the function is of type 'void', it:
        1) replaces each return statement inside the new sequence with
           'goto __exit',
        2) appends an empty exiting label '__exit' to the sequence;
    if the function returns a value, it:
        1) preprends the declaration of a return variable '__res' to the
           sequence,
        2) replaces each return statement inside the sequence with
           '__res = x; goto __exit'.
        3) appends the final and unique labelled return statement
           '__exit; return __res;' to the sequence.

    [mark] - mark to put on the sequence the function's body is wrapped into,
    [tg] - target function definition AST term. *)
(*========================*)
(** [unfold_let_mult_on ?constify t]: see [unfold_let_mult]. *)
(*========================*)
(** [unfold_let_mult ~constify tg]: expects target [tg] to point
   at a multiple variable declaration. Then, it replaces it by a sequence of
   simple variable declarations.

   For example:

     int a = 1, b = a;

   becomes:

     int a = 1;
     int b = a;

   If [constify] is [true] one or more variable declarations from the target
   multiple variable declaration [tg] will be constified. See comments in
   [unfold_let_mult_on] for more details. *)
(*========================*)
(** [constify_args_on ?force t]: see [constify_args]. *)
(*========================*)
(** [constify_args ?force tg]: expects the target [tg] to point at a function
   definition. Then, based on the constification records in [const_records], it
   shall constify the arguments that should be constified. If the corresponding
   constification record says so, it shall also constify the target function
   itself.

   One can force, e.g. for testing purposes, the function to ignore the
   constification records and to constify all of the function's arguments as
   well as the function itself. *)
(*========================*)
(** [constify_aliases_on ?force t]: see [constify_aliases]. *)
(*========================*)
(** [constify_aliases ?force tg]: expects target the target [tg] to point at a
   function definition. Then, based on the constification records in
   [const_records], it shall constify the aliases to arguments that have been
   constified during the constification process.

   One can force, e.g. for testing purposes, the function to ignore the
   constification records and to constify the aliases to all of the function's
   arguments. *)
(*========================*)
(** [heapify_on t]: see [heapify]. *)
(*========================*)
(** [heapify tg]: expects the target [tg] to point at a simple or a multiple
   variable declaration. Then, if it is necessary, i.e. if the variable is not a
   reference or a pointer to a previously user-allocated memory location, the
   function shall promote the variable from the stack to the heap.

   Example:

   int tab[5] = { 1, 2, 3, 4, 5 };

   becomes:

   int * tab = new int[5] { 1, 2, 3, 4, 5 }

   However, in:

   int * a = new int(10);
   int &b = &a;

   nor [a] nor [b] are transformed. *)
(*========================*)
(** [vars_tbl]: hashtable generic to keep track of variables and their pointer
    depth. This abstrastion is used for generic functions. *)
(*========================*)
(** [get_vars_data_from_cptr_arith va t] : resolve pointer operation to get the
    pointer variable. Then, return the data of the corresponding variable stored
    in vars_tbl. *)

(* //////////////////////////////
   //////  Align_basic.ml  //////
   ////////////////////////////// *)

      (*========================*)
(** [def vec_align tg]: expects the target [tg] to point at an array declaration, then it will add the alignas attribute
    to its type with [vec_align] alignment size. *)
(*========================*)
(** [header ()]: inserts "#include \"stdalign.h\"" at the top of the main file. *)
(*========================*)
(** [assume alignment var tg]: expects the target [tg] to be a relative target, then it will insert an instruction of the
    form var = __builtin_assume_aligned(var, alignment) at that location.*)

(* //////////////////////////////
   /////  Rewrite_core.ml  //////
   ////////////////////////////// *)

      (*========================*)
(** [apply_rule_on rule t]: applies rule [rule] on trm [t]. *)
(*========================*)
(** [compute_on t]: applies arithmetic simplifications on trm [t]. *)

(* //////////////////////////////
   /////////  Arith.ml  /////////
   ////////////////////////////// *)

      (*========================*)
(** [simpl_surrounding_expr] first goes to the outside of the targeted expression,
   then applies [simpl] *)

(* //////////////////////////////
   /////  Typedef_core.ml  //////
   ////////////////////////////// *)

      (*========================*)
(** [fold_at fold_at index]: replaces occurrences of the typedef underlying type with the defined type,
       [fold_at] - targets where folding should be performed, if left empty then folding is applied recursively
                 on all the trms that belong to the same sequence as typedef,
       [index] - index of the typedef on its surrounding sequence,
       [t] - ast of the sequence that contains the targeted typedef. *)
(*========================*)
(** [unfold_at unfold_at]: replaces occurrences of the defined type with its underlying type,
      [delete] - a flag for deciding if we should delete or not the typedef declaration,
      [unfold_at] - targets where unfolding should be performed, if empty unfolding is applied
                   on all the trms that are with in the same level as the targeted typedef,
      [t] - ast of the sequence that contains the targeted typedef. *)
(*========================*)
(** [insert_copy_of name index t]: creates a copy of the targeted typedef
      [name] - new typ name
      [t] - ast of the surrounding sequence of the targeted typedef *)
(*========================*)
(** [insert_at name td_body index]: inserts a new type definition,
      [name] - new type name,
      [td_body] - body of the new type definition,
      [index] - location where the typedef should be inserted inside a sequence. *)

(* //////////////////////////////
   //////  Cast_basic.ml  ///////
   ////////////////////////////// *)

      (*========================*)
(** [insert ty tg]: expects the target [tg] to point at any expression that can be casted,
    then it will cast it to type [ty] *)

(* //////////////////////////////
   //////  Instr_basic.ml  //////
   ////////////////////////////// *)

      (*========================*)
(** [delete tg]: expects the target [tg] to point at an instruction inside a sequence
      then it will remove that instruciton from that sequence *)
(*========================*)
(** [copy ~target tg]: expects the target [tg] to point at an instruction that is
    going to be copied to the target [dest].
    TODO: Either support copying at arbitrary location or use a relative target *)
(*========================*)
(** [move ~rev ~dest tg] move the instructions at target [tg] to the relative position [dest].
    If you want to move contiguous instructions together, you can use a span target.

    The relative position [dest] is computed independently for each target so it is
    allowed to have matches inside the same sequence, inside different sequences,
    or a mix between the two as long as you can write a common [dest] target.

    If [rev] is true, perform the moves in the reverse order if the target resolves
    to multiple paths.

    TODO: decide and specify how marks between instructions are moved in the sequence

   @correctness: Correct if the swapped instructions are parallelizable:
   {H1 * H} instr1 {H1' * H} and {H2 * H} instr2 {H2' * H}
   which lead globally to the derivation
   {H1 * H2 * H} instr1 {H1' * H2 * H} instr2 {H1' * H2' * H}
   we can build the same postcondition with
   {H1 * H2 * H} instr2 {H1 * H2' * H} instr1 {H1' * H2' * H}

   This is sufficient but not necessary, a manual commutation proof can be used
   as well. *)
(*========================*)
(** [read_last_write ~write tg]: expects the target [tg] to point at a read operation,
    then it replaces the read operation with left hand side of the write operation targeted
    by [write]

   @correctness: the read expression must be pure, and its evaluation must not
   have changed since the write.*)
(*========================*)
(** [accumulate tg] expects the target [tg] to point at a block of write operations that write to the same memory location
    then accumulate those write operations into a single one.
    Ex.
    int x;
    {
      x += 1;
      x += 2;
      x += 3;
    }
    is transformed to x += 1+2+3 *)

(* //////////////////////////////
   ///////  Loop_core.ml  ///////
   ////////////////////////////// *)

      (*========================*)
(** [color_on nb_colors i_color t]: transform a loop into two nested loops based on the coloring pattern,
      [nb_colors] - a variable used to represent the number of colors,
      [i_color] - a variable representing the index used of the new outer loop,
      [t] - ast of the loop. *)
(*========================*)
(** [tile_bound]: used for loop tiling transformation *)
(*========================*)
(** [tile_on divides b tile_index t]: tiles loop [t],
      [tile_index] - string representing the index used for the new outer loop,
      [bound] - a tile_bound type variable representing the type of the bound used in
                 this transformation,
      [t] - ast of targeted loop. *)
(*========================*)
(** [fusion_on_block_on t]: merges two or more loops with the same components except the body,
      [t] - ast of the sequence containing the loops. *)
(*========================*)
(** [grid_enumerate_on indices_and_bounds t]: transforms a loop over a grid into nested loops over
    each dimension of that grid,
      [indices_and_bounds] - a list of pairs representing the index and the bound for each dimension,
      [t] - ast of the loop. *)
(*========================*)
(** [unroll_on index t]: unrolls loop [t],
      [inner_braces] - a flag on the visibility of generated inner sequences,
      [outer_seq_mark] - generates an outer sequence with a mark,
      [t] - ast of the loop. *)
(*========================*)
(** [unswitch_at trm_index t]: extracts an if statement inside the loop whose condition,
    is not dependent on the index of the loop or any other local variables,
      [trm_index] - index of the if statement inside the body of the loop,
      [t] - ast of the for loop. *)
(*========================*)
(** [to_unit_steps_on new_index t]: transforms loop [t] into a loop with unit steps,
      [new_index] - a string representing the new index for the transformed loop,
      [t] - ast of the loop to be transformed. *)
(*========================*)
(** [scale_range factor new_index t]: transforms loop [t] into a loop with larger steps,
      [factor] - the multiplicative factor on the indices,
      [new_index] - a string representing the new index for the transformed loop,
        (if empty, then we reuse the original name)
      [t] - ast of the loop to be transformed. *)
(*========================*)
(** [fold_at index start step t]: transforms a sequence of instructions into a for loop,
      [index] - index of the generated for loop,
      [start] - starting value for the index of the generated for loop,
      [step] - step of the generated for loop,
      [t] - ast of the sequence.

    NOTE: we trust the user that "stop" corresponds to the number of iterations
    LATER: use  sExpr  to mark the subexpression that correspnod to the string "start";
    then you can Generic.replace at these marks.*)
(*========================*)
(** [split_range_at nb cut]: splits a loop into two loops based on the range,
     [nb] - by default this argument has value 0, if provided it means that it will split the loop at start + nb iteration,
     [cut] - by default this argument has value tmr_unit(), if provided then the loop will be splited at that iteration,
     [t] - ast of the for loop.

     TODO: Optional arguments instead of weird "default" values
     *)
(*========================*)
(** [rename_index_on new_index]: renames the loop index variable *)

(* //////////////////////////////
   ///////  Resources.ml  ///////
   ////////////////////////////// *)

      (*========================*)
(** [trm_is_pure]: Does this term always evaluate to the same value?
    The computed value should no be affected by observable program state.
    Computing the value should not affect the observable program state. *)
(*========================*)
(** [fun_minimize]: minimize a function contract by looking at the resource usage of its body *)
(*========================*)
(** [loop_minimize]: minimize linear invariants of a loop contract *)
(*========================*)
(** [make_strict_loop_contracts] uses computed resources to fix the default loop contracts on all the children of the targeted node *)
(*========================*)
(** [detach_loop_ro_focus tg] transforms all the ressources that are in a reads clause into a resource in a par_read clause with a focus around the loop body. *)
(*========================*)
(** [specialize_arbitrary_fracs tg] specializes all the arbitrarily chosen fractions that appear at the given point in a sequence.

  If a given ressource has n carved holes at the split point, then each of these will be specialized to 1/(n+1) and ghosts will be inserted to pre-split the resource as early as possible.

  Currently this does not check that the eliminated arbitrary fractions are not used outside the sequence immediately surrounding the target. This may be a cause of typing errors after calling this transformation.
*)

(* //////////////////////////////
   /////////  Expr.ml  //////////
   ////////////////////////////// *)

      (*========================*)
(** [replace_fun code tg]: expects the target [tg] to point at a function call,
    then it replaces the name of the function call with the one entered by the user

    Assumption:
      [name] is the name of an already defined function which has the same signature
      as function whose call is targeted by [tg] *)

(* //////////////////////////////
   //////  Matrix_core.ml  //////
   ////////////////////////////// *)

      (*========================*)
(** [tile_none]: special tile value where the dimension should be kept fully *)
(*========================*)
(** [tile_none]: special tile value where the dimension should be dropped *)
(*========================*)
(** [memcpy dest d_offset d_dims src s_offset s_dims copy_dims elem_size] uses a series
    of [matrixN_contiguous] and [mindexN_contiguous] ghosts to issue a [MMEMCPY] from
    matrix [dest] at [d_offset] to matrix [src] at [s_offset].

    Preconditions:
    - [dest] has dims [d_dims] and [(len d_offset) <= (len d_dims)]
    - [src] has dims [s_dims] and [(len s_offset) <= (len s_dims)]
    - [(len d_dims) - (len d_offset) = (len s_dims) - (len s_offset) = (len copy_dims)]
    - [take_last (len copy_dims) d_dims = take_last (len copy_dims) s_dims]
  *)
(*========================*)
(** [map_all_accesses v dims map_indices mark t] maps all accesses to [v] in [t],
    using the [map_access] function.

    Fails if [v] occurs in a sub-term that is not an access to [v], as this could mean some
    accesses are hidden (e.g. behind a function call).

    - the dimensions and type of the matrix [v] are stored in [ret_dims_and_typ] if provided.
    *)
(*========================*)
(** [replace_all_accesses prev_v v dims map_indices mark t] replaces all accesses to [prev_v]
    in [t] with accesses to [v], using new [dims] and changing indices with [map_indices].
    *)
(*========================*)
(** [pointwise_fors ?reads ?writes ?modifies ranges body] creates nested loops
  with [ranges] over the main body [body].
  The body has the given [reads], [writes], and [modifies].
  Each loop contract adds a layer of pointwise Group resources.
  *)
(*========================*)
(** [intro_calloc_aux t]: replaces a call to calloc with a call to the macro CALLOC,
     [t] - ast of the call to alloc. *)
(*========================*)
(** [intro_malloc_aux t]: replaces a call to calloc with a call to MALLOC,
     [t] - ast of the call to alloc. *)
(*========================*)
(** [intro_mindex_aux dim t] replaces an array access at index [i] with an array access at MINDEX([dim],i),
      [dim] - the size of the array accesses with [t],
      [t] - the ast of the array access. *)
(*========================*)
(** [reorder_dims_aux order t]: reorders the dimensions in a call to CALLOC, MALLOC or MINDEX,
      [order] - a list of indices based on which the elements in dims should be ordered,
      [t] - ast of the call to CALLOC, MALLOC, MINDEX. *)
(*========================*)
(** [insert_alloc_dim_aux new_dim t]: adds a new dimension at the beginning of the list of dimension,
     [new_dim]: the new dimension which is goin to be inserted into the list of dims in call to CALLOC or MALLOC,
     [t]: ast of the call to ALLOC functions. *)
(*========================*)
(** [insert_access_dim_index_aux new_dim new_index t]: add a new dimension at the beginning of the list of dimension
     and add a new index at the begining of the list of indices in the call to MINDEX inside the
     targeted array access.
      [new_dim]: the new dimension which is goin to be inserted into the list of dims in call to MINDEX,
      [new_index]: the new index which is goin to be inserted into the list of indices in call to MINDEX,
      [t]: ast of the array_access with the updated list of args in the call to MINDEX. *)
(*========================*)
(** [local_name_aux mark var local_var malloc_trms var_type t] insert a local matrix declaration with name [local_var] and copy the content
      from the matrix [var] to [local_var] and replace all the accesses of that matrix inside the targeted insturction [t].
      [mark] - an optional mark at the final generated sequence,
      [var] - the name of the current matrix used in instruction [t],
      [new_var] - the name of the local matrix which replaces all the current accesses of [var],
      [t] - ast of thee instuction which contains accesses to [var]. *)

(* //////////////////////////////
   /////  Variable_core.ml  /////
   ////////////////////////////// *)

      (*========================*)
(** [fold_decl_at fold_at index t]: fold the targeted variable definition,
      [fold_at] - target where folding should be performed, if left empty
                  then folding is applied everywhere,
      [index] - the index of the targeted definition on its surroudinig sequence,
      [t] - ast of the sequence that contains the targeted definition. *)
(*========================*)
(** [rename_at new_name index t]: renames all the occurences of the variable declared on the targeted declaration,
      [new_name] - the new name for the targeted variable,
      [index] - index of the targeted declaration inside its surrounding sequence,
      [t] - ast of the sequence that contains the targeted declaration. *)
(*========================*)
(** [init_detach_on t]: detaches the targeted variable declaration,
      [t] - ast of the targeted variable declaration. *)
(*========================*)
(** [Init_attach_no_occurrences]: raised by [init_attach_at]. *)
(*========================*)
(** [Init_attach_occurrence_below_control]: raised by [init_attach_at]. *)
(*========================*)
(** [init_attach_at t]: attaches a variable declaration to its unique write operation,
      [const] - a boolean to decide if the attached variable should be mutable or not,
      [index] - index of the targeted instruction inside its surrounding sequence,
      [t] - ast of the surrounding sequence of the variable declaration.

    NOTE: if no set operation on the targeted variable was found then Init_attach_no_occurrences is raised
          if more then one set operation on the targeted variable was found then Init_attach_occurrence_below_control is raised *)
(*========================*)
(** [delocalize_at array_size ops index t]: see [Variable_basic.delocalize],
      [array_size] - size of the arrays to be declared inside the targeted sequence,
      [ops] - delocalize operation representing the unitary lement used for initialization
             and the fold_lefting operation used for the reduction,
      [index] - the index for the two added loops,
      [t] - the ast of the sequence generated after applying the local_name transformation. *)
(*========================*)
(** [insert_at index const name typ value t]: inserts a variable declaration on sequence [t],
      [index] - location where the declaration is going to be inserted,
      [const] - a flag on the mutability of the variable [name],
      [name] - name of the inserted variable,
      [typ] - the type of the inserted variable,
      [value] - the initial value of the inserted variable [name] entered as a string,
      [t] - ast of the sequence where the insertion is performed. *)
(*========================*)
(** [change_type_at new_type t]: changes the current type of the targeted variable,
      [new_type] - the new type replacing the current one entered as a string,
      [t] - ast of the sequence that contains the targeted declaration. *)
(*========================*)
(** [bind_at index fresh_name const p_local t]: binds the variable [fresh_name] to the targeted trm,
      [mark_let] - an optional mark attached to the new binding instruction
      [mark_occ] - an optional mark attached to the occurrences of [fresh_name] that are inserted
      [my_mark] - a mark to be left on the bound term, inside the new let-binding definition
      [index] - index of the instruction containing the targeted function call,
      [fresh_name] - name of the variable which going to be binded to the function call,
      [const] - a flag for the mutability of the binded variable,
      [p_local] - the local path from the instruction containing the targeted node
                  to the targeted node,
      [t] - ast of the sequence containing the targeted node. *)
(*========================*)
(** [remove_get_operations_on_var x t]: removes one layer of get operations on variable [x].
     i.e. if [x] was a pointer to [v], [get x] becomes [v].
   *)
(*========================*)
(** [remove_get_operations_on_var_temporary x t]: to be removed. *)
(*========================*)
(** [to_nonconst_at index t]: transforms a constant into a mutable variable.
      [index] - the index of the targeted declaration inside its surrounding sequence,
      [t] - ast of the sequence that contains the targeted declaration. *)
(*========================*)
(** [to_const_at index t]: transform a mutable variable without explicit writes into a constant,
      [index] - the index of the targeted declaration inside its surrounding sequence,
      [t] - ast of the sequence that contains the targeted declaration. *)
(*========================*)
(** [simpl_deref_on t]: checks if [t] is of the form *(&b) or *(&b),
    if that's the case then simplify that expression and return it,
       [indepth] - search indepth for the targeted expressions,
       [t] - trm that represents one of the epxressions *(&b) or &( *b). *)
(*========================*)
(** [ref_to_pointer_at index t]: transforms the targeted declaration from a reference to a poitner,
      [index] - index of that targeted declaration in its surrounding block,
      [t] - ast of the sequence that contains the targeted declaration. *)
(*========================*)
(** [ref_to_var_on t]: converts a reference variable to a simple stack var variable
     [t] - ast of the refernce declaration *)

(* //////////////////////////////
   /////  Accesses_core.ml  /////
   ////////////////////////////// *)

      (*========================*)
(** [transform_on f_get f_set t]: applies f_get or f_set depending on the fact if
    [t] is a get operation or a set operation,
      [f_get] - the get operation that is going to be applied,
      [f_set] - the set operation that is going to be applied,
      [t] - the ast of the node where the operation is applied to. *)
(*========================*)
(** [intro_on t]: changes encodings "struct_get(get (t), f)" to "get(struct_access (t, f))",
      [t] - ast of the node where the accesses can be found. *)

(* //////////////////////////////
   ////////  Matrix.ml  /////////
   ////////////////////////////// *)

      (*========================*)
(** [intro_calloc tg]: expects the target [tg] to point at a variable declaration
    then it will check its body for a call to calloc. On this extended path it will call
    the [Matrix_basic.intro_calloc] transformation *)
(*========================*)
(** [intro_mindex dim]; expects the target [tg] to point at at a matrix declaration, then it will change
     all its occurrence accesses into Optitrust MINDEX accesses. *)
(*========================*)
(** [intro_malloc tg]: expects the target [tg] to point at a variable declaration
    then it will check its body for a call to malloc. On this extended path it will call
    the [Matrix_basic.intro_malloc] transformation. *)
(*========================*)
(** [biject fun_bij tg]: expects the target [tg] to point at at a matrix declaration , then it will search for all its
    acccesses and replace MINDEX with  [fun_bij]. *)
(*========================*)
(** [intro_mops dims]: expects the target [tg] to point at an array declaration allocated with
      calloc or malloc, then it will apply intro_calloc or intor_mmaloc based on the type of
      the current allocation used. Then it will search for all accesses and apply intro_mindex. *)
(*========================*)
(** [elim_mops]: expects the target [tg] to point at a subterm and
  eliminates all MINDEX macros in that subterm.

  TODO:
  - eliminate MALLOC2 into malloc(sizeof(T[n][m]))?
  - ~simpl
*)
(*========================*)
(** [delocalize ~mark ~init_zero ~acc_in_place ~acc ~last ~var ~into ~dim ~index ~indices ~ops tg]: this is a combi
   varsion of [Matrix_basic.delocalize], this transformation first calls Matrix_basi.local_name to create the isolated
    environment where the delocalizing transformatino is going to be performed *)
(*========================*)
(** [reorder_dims ~rotate_n ~order tg] expects the target [tg] to point at at a matrix declaration, then it will find the occurrences of ALLOC and INDEX functions
      and apply the reordering of the dimensions. *)
(*========================*)
(** [elim]: eliminates the matrix [var] defined in at the declaration targeted by [tg].
  All reads from [var] must be eliminated The values of [var] must only be read locally, i.e. directly after being written.
  *)
(*========================*)
(** [inline_constant]: expects [tg] to target a matrix definition,
   then first uses [Matrix.elim_mops] on all reads before attempting
   to use [Arrays.inline_constant].
   *)
(*========================*)
(** [elim_constant]: expects [tg] to target a matrix definition,
   then first uses [Matrix.elim_mops] on all reads before attempting
   to use [Arrays.elim_constant].
   *)
(*========================*)
(** [iter_on_var_defs]: helper for transformations that need to iterate
  on variable definitions while requiring the path to the surrounding sequence.
   *)
(*========================*)
(** [delete] expects target [tg] to point to a definition of matrix [var], and deletes it.
  Both allocation and de-allocation instructions are deleted.
  Checks that [var] is not used anywhere in the visible scope.
   *)
(*========================*)
(** [local_name_tile ?delete ?indices ?alloc_instr ?local_var tile ?simpl tg] is a convenient
  version of {!Matrix_basic.local_name_tile}. It deletes the original matrix if [delete = true]
  or if [local_var = ""].
   *)

(* //////////////////////////////
   /////  Arrays_basic.ml  //////
   ////////////////////////////// *)

      (*========================*)
(** [to_variables new_vars tg]: expects the target [tg] to point at an array declaration.
    Then it transforms this declaration into a list of declarations.
    [new_vars] - denotes the list of variables that is going to replace the initial declaration
      the length of this list is equal to [size -1] where [size] is the size of the array.*)
(*========================*)
(** [tile ~block_type block_size tg]: expects the target [tg] to point at an array declaration.
   Then it takes that declaration and transforms it into a tiled array. All the accesses of the
   targeted array are handled as well.
   [block_type] - denotes the name of the array which is going to represent a tile.
   [block_size] - size of the block of tiles. *)
(*========================*)
(** [swap name x tg]: expects the target [tg] to point at an array declaration.
   It changes the declaration so that the bounds of the array are switched. Also
   all the accesses of the targeted array are handled as well.*)
(*========================*)
(** [aos_to_soa tv sz] finds the definition of type [tv] which should be a typedef Record.
    Then it will change its struct fields type to arrys of size [sz] with type their current type.
    All the accesses will be swapped.
    Ex:
      int const N = 100;
      typedef struct {
        int x;
        int y;
      } vect;
      vect w[N];
      int main(){
        int i;
        int c = w[i].x;
        return 0;
      }

      int const N = 100;
      typedef struct {
        int x[N];
        int y[N];
      }
      vect w
      int main(){
        int i;
        int c = w.x[i];
        return 0;
      }
*)
(*========================*)
(** [set_explicit tg] expects the target [tg] to point at an array declaration
    then it will remove the initialization trm and a list of write operations on
    each of the cells of the targeted array.
*)
(*========================*)
(** [inline_constant] expects the target [decl] to point at a constant array literal declaration, and resolves all accesses targeted by [tg], that must be at constant indices.
  *)
(*========================*)
(** [elim] expects the target [tg] to point at a constant array literal declaration, and eliminates it if it is not accessed anymore.
  *)

(* //////////////////////////////
   ///  Specialize_basic.ml  ////
   ////////////////////////////// *)

      (*========================*)
(** [any e tg]: expects the target [tg] to be point at a call to the function [ANY], then it will replace it with [e]. *)
(*========================*)
(** [choose_fct select_arg]: expects the target [tg] to point at a call to the function [CHOOSE]
    which is used in the delocalize transformation (see to [Variable.delocalize] ),
    then it will replace that call with one of its arguments that satisfies the predicate [select_arg]. *)
(*========================*)
(** [choose_id id tg]: chooses the id of the arguments of the function [CHOOSE], then this id is used
    by the function [choose_fct]. *)
(*========================*)
(** [choose choice tg]: combines [choose_fct] and [choose_id] into one function so that [choice] is used
    when applying function [choose_fct]. *)
(*========================*)
(** [fundefs spec_name spec_args tg] *)
(*========================*)
(** [funcalls spec_name args_to_choose tg]: expects the target [ŧg] to point to a function call, and assumes that
      there is already a function generated by the transformation [fundefs], then it will replace that
      call with a call to the function [spec_name] already defined either by using the transformation [fundefs],
      Or manually by the user.*)

(* //////////////////////////////
   //////  Instr_core.ml  ///////
   ////////////////////////////// *)

      (*========================*)
(** [copy_at dest_index index t]: copies instruction at [index] to the [dest_index],
     [dest_index] - the [index] where the target instruction will be copied,
     [index] - index of the targeted instruction,
     [t] - ast of the surrounding sequence of the targeted instruction. *)
(*========================*)
(** [move_at dest_index index t]: moves instruction at [index] to the [dest_index],
     [dest_index] - the [index] where the target instruction will be copied,
     [index] - index of the targeted instruction,
     [t] - ast of the surrounding sequence of the targeted instruction. *)
(*========================*)
(** [accumulate_on t]: transform a list of write instructions into a single instruction,
      [t] - the ast of the sequence containing the instructions. *)

(* //////////////////////////////
   //////  Expr_basic.ml  ///////
   ////////////////////////////// *)

      (*========================*)
(** [update f tg]: applies the operation [f] at the target [tg] *)
(*========================*)
(** [replace node tg]: expects the target to point at an instruction, then it will replace this
    instruction with [node]. Note that [node] can be also some code entered as string if that is
    the case then to integrate it on the current ast this transformation shoudl be called with the flag ~reparse:true

   @correctness: Needs local manual reproving that if an invariant in the
   previous proof was { H } old_expr { Q } then { H } new_expr { Q } holds
   as well *)
(*========================*)
(** [replace_fun code tg]: expects the target to point at a function call,
    then it replaces the name of the function call with the one entered
    by the user

    Assumption:
      [name] is the name of an already defined function which has the same
      signature as function whose call is targeted by [tg] *)
(*========================*)
(** [view_subterms tg]: displays on stdout all the subterms of the targeted term.
   For viewing on stdout all subterms of a program, use:
     Expr.view_subterms [];
   which is like
     Expr.view_subterms [dRoot];
   and possibly specify a regexp to investigate:
     Expr.view_subterms ~constr:(sInstr "+= 2") [dRoot];
   Note that this is for debugging purpose only. *)
