(**

This file describes transformations of the layout of records and arrays.

Author: Ramon Fernandez I Mir and Arthur Charguéraud.

License: MIT.

*)

Set Implicit Arguments.
Require Export Semantics LibSet LibMap TLCbuffer Typing.



(* ********************************************************************** *)
(* * Specification of the transformation *)

Module Example.

(* Initial typdefctx *)
Definition pos : map field typ := (\{})["x" := typ_int]["y" := typ_int]["z" := typ_int].
Definition C : typdefctx := (\{})["pos" := (typ_struct pos)].

(* Final typdefctx *)
Definition struct_x : map field typ := (\{})["x" := typ_int].
Definition pos' : map field typ := (\{})["s" := (typ_var "struct_x")]["y" := typ_int]["z" := typ_int].
Definition C' : typdefctx := (\{})["pos" := (typ_struct pos')]["struct_x" := (typ_struct struct_x)].

End Example.


(* ********************************************************************** *)
(* * Definition of the transformation *)

(** Grouping transformation. Specified by:
    - The name of the struct to be modified.
    - The set of fields to be grouped.
    - The name of the new struct that will hold the fields to be grouped.
    - The name of the field in the new struct that will have as type
      the new struct. *)

Record group_tr := make_group_tr {
  group_tr_struct_name : typvar;
  group_tr_fields : set field;
  group_tr_new_struct_name : typvar;
  group_tr_new_struct_field : field
}.

Notation make_group_tr' := make_group_tr.

(** Checking if the transformation is acceptable *)

Inductive group_tr_ok : group_tr -> typdefctx -> Prop :=
  | group_tr_ok_intros : forall Tfs Tt fs fg Tg gt C,
      gt = make_group_tr Tt fs Tg fg ->
      Tt \indom C ->
      C[Tt] = typ_struct Tfs ->
      Tg \notindom C ->
      fs \c dom Tfs ->
      fg \notindom Tfs ->
      group_tr_ok gt C.


(* ********************************************************************** *)
(* * The transformation applied to the different constructs. *)

(** Transformation of typdefctxs: C ~ |C| *)

Inductive tr_typdefctx (gt:group_tr) : typdefctx -> typdefctx -> Prop :=
  | tr_typdefctx_intro : forall Tfs Tfs' Tfs'' Tt fs Tg fg C C',
      gt = make_group_tr Tt fs Tg fg ->
      dom C' = dom C \u \{Tg} ->
      C[Tt] = typ_struct Tfs ->
      C'[Tt] = typ_struct Tfs' ->
      C'[Tg] = typ_struct Tfs'' ->
      (forall T,
        T \indom C ->
        T <> Tt ->
        C'[T] = C[T]) ->
      dom Tfs' = (dom Tfs \- fs) \u \{fg} ->
      (forall f,
        f \indom Tfs ->
        f \notin fs ->
        Tfs'[f] = Tfs[f]) ->
      dom Tfs'' = fs ->
      (forall f,
        f \indom Tfs'' ->
        Tfs''[f] = Tfs[f]) ->
      tr_typdefctx gt C C'.

(** Transformation of paths: π ~ |π| *)

Inductive tr_accesses (gt:group_tr) : accesses -> accesses -> Prop :=
  | tr_accesses_nil :
      tr_accesses gt nil nil
  | tr_accesses_array : forall π π' T i,
      tr_accesses gt π π' ->
      tr_accesses gt ((access_array T i)::π) ((access_array T i)::π')
  | tr_accesses_field_group : forall Tt fs fg Tg f a0 π a1 a2 π',
      tr_accesses gt π π' ->
      gt = make_group_tr Tt fs Tg fg ->
      f \in fs ->
      a0 = access_field (typ_var Tt) f ->
      a1 = access_field (typ_var Tt) fg ->
      a2 = access_field (typ_var Tg) f ->
      tr_accesses gt (a0::π) (a1::a2::π')
  | tr_accesses_field_other : forall T Tt π π' f,
      tr_accesses gt π π' ->
      Tt = group_tr_struct_name gt ->
      (T <> (typ_var Tt) \/ f \notin (group_tr_fields gt)) ->
      tr_accesses gt ((access_field T f)::π) ((access_field T f)::π').

(** Transformation of values: v ~ |v| *)

Inductive tr_val (gt:group_tr) : val -> val -> Prop :=
  | tr_val_uninitialized :
      tr_val gt val_uninitialized val_uninitialized
  | tr_val_unit :
      tr_val gt val_unit val_unit
  | tr_val_bool : forall b,
      tr_val gt (val_bool b) (val_bool b)
  | tr_val_int : forall i,
      tr_val gt (val_int i) (val_int i)
  | tr_val_double : forall d,
      tr_val gt (val_double d) (val_double d)
  | tr_val_abstract_ptr : forall l π π',
      tr_accesses gt π π' ->
      tr_val gt (val_abstract_ptr l π) (val_abstract_ptr l π')
  | tr_val_array : forall a T a',
      length a = length a' ->
      (forall i,
        index a i ->
        tr_val gt a[i] a'[i]) ->
      tr_val gt (val_array T a) (val_array T a')
  | tr_val_struct_group : forall Tt Tg s s' fg fs sg,
      gt = make_group_tr Tt fs Tg fg ->
      fs \c dom s ->
      fg \notindom s ->
      dom s' = (dom s \- fs) \u \{fg} ->
      dom sg = fs ->
      (forall f,
        f \indom sg ->
        tr_val gt s[f] sg[f]) ->
      (forall f,
        f \notin fs ->
        f \indom s ->
        tr_val gt s[f] s'[f]) ->
      s'[fg] = val_struct (typ_var Tg) sg ->
      tr_val gt (val_struct (typ_var Tt) s) (val_struct (typ_var Tt) s')
  | tr_val_struct_other : forall Tt T s s',
      Tt = group_tr_struct_name gt ->
      T <> (typ_var Tt) ->
      dom s = dom s' ->
      (forall f,
        f \indom s ->
        tr_val gt s[f] s'[f]) ->
      tr_val gt (val_struct T s) (val_struct T s').

Definition is_struct_op (op:prim) :=
  match op with
  | prim_struct_access _ _ => True
  | prim_struct_get _ _ => True
  | _ => False
  end.

(* Transformation used in the struct cases to avoid repetition. *)

Inductive tr_struct_op (gt:group_tr) : trm -> trm -> Prop :=
  | tr_struct_op_group : forall fs Tt fg pr Tg f op1 op2 ts,
      pr = prim_struct_access \/ pr = prim_struct_get ->
      gt = make_group_tr Tt fs Tg fg ->
      f \in fs ->
      op1 = pr (typ_var Tt) fg ->
      op2 = pr (typ_var Tg) f ->
      tr_struct_op gt (trm_app (pr (typ_var Tt) f) ts) (trm_app op2 ((trm_app op1 ts)::nil))
  | tr_struct_op_other : forall Tt T pr f ts,
      pr = prim_struct_access \/ pr = prim_struct_get ->
      Tt = group_tr_struct_name gt ->
      (T <> (typ_var Tt) \/ f \notin (group_tr_fields gt)) ->
      tr_struct_op gt (trm_app (pr T f) ts) (trm_app (pr T f) ts).

(** Transformation of terms: t ~ |t| *)

Inductive tr_trm (gt:group_tr) : trm -> trm -> Prop :=
  | tr_trm_val : forall v v',
      tr_val gt v v' ->
      tr_trm gt (trm_val v) (trm_val v')
  | tr_trm_var : forall x,
      tr_trm gt (trm_var x) (trm_var x)
  | tr_trm_if : forall t1 t2 t3 t1' t2' t3',
      tr_trm gt t1 t1' ->
      tr_trm gt t2 t2' ->
      tr_trm gt t3 t3' ->
      tr_trm gt (trm_if t1 t2 t3) (trm_if t1' t2' t3')
  | tr_trm_let : forall x t1 t2 t1' t2',
      tr_trm gt t1 t1' ->
      tr_trm gt t2 t2' ->
      tr_trm gt (trm_let x t1 t2) (trm_let x t1' t2')
  (* Special case: structs *)
  | tr_trm_struct_op : forall t1' op t1 tr,
      is_struct_op op ->
      tr_trm gt t1 t1' ->
      tr_struct_op gt (trm_app op (t1'::nil)) tr ->
      tr_trm gt (trm_app op (t1::nil)) tr
  (* Args *)
  | tr_trm_args0 : forall op,
      tr_trm gt (trm_app op nil) (trm_app op nil)
  | tr_trm_args1 : forall op t1 t1',
      ~ is_struct_op op ->
      tr_trm gt t1 t1' ->
      tr_trm gt (trm_app op (t1::nil)) (trm_app op (t1'::nil))
  | tr_trm_args2 : forall op t1 t1' t2 t2',
      tr_trm gt t1 t1' ->
      tr_trm gt t2 t2' ->
      tr_trm gt (trm_app op (t1::t2::nil)) (trm_app op (t1'::t2'::nil)).

(** Transformation of stacks: S ~ |S| *)

Inductive tr_stack_item (gt:group_tr) : (var * val) -> (var * val) -> Prop :=
  | tr_stack_item_intro : forall x v v',
      tr_val gt v v' -> 
      tr_stack_item gt (x, v) (x, v').

Inductive tr_stack (gt:group_tr) : stack -> stack -> Prop :=
  | tr_stack_intro : forall S S',
      LibList.Forall2 (tr_stack_item gt) S S' ->
      tr_stack gt S S'.

Lemma stack_lookup_tr : forall gt S S' x v,
  tr_stack gt S S' ->
  Ctx.lookup x S = Some v -> 
    exists v', 
       Ctx.lookup x S' = Some v' 
    /\ tr_val gt v v'.
Proof.
  introv HS Hx. inverts HS as HS. induction HS.
  { inverts Hx. }
  { inverts H as Hv. inverts Hx as Hx. case_if in Hx.
    { inverts Hx. exists v'. splits*. unfolds. case_if*. }
    { forwards (v''&Hx'&Hv''): IHHS Hx. exists v''.
      splits*. unfolds. case_if. fold Ctx.lookup. auto. } }
Qed.

(** Transformation of states: m ~ |m| *)

Inductive tr_state (gt:group_tr) : state -> state -> Prop :=
  | tr_state_intro : forall m m',
      dom m = dom m' ->
      (forall l,
        l \indom m ->
        tr_val gt m[l] m'[l]) ->
      tr_state gt m m'.


(* ---------------------------------------------------------------------- *)
(** Hints *)

(* TODO: Factor this out. *)
Lemma index_of_index_length : forall A (l' l : list A) i,
  index l' i ->
  length l' = length l ->
  index l i.
Proof.
  intros. rewrite index_eq_index_length in *.
  applys* index_of_index_length'.
Qed.

Hint Resolve index_of_index_length.

Hint Constructors red.
Hint Constructors tr_trm tr_val tr_accesses tr_state tr_stack.
Hint Constructors read_accesses write_accesses.


(* ---------------------------------------------------------------------- *)
(** Functionality of the relations *)

Theorem functional_tr_accesses : forall gt π π1 π2,
  tr_accesses gt π π1 ->
  tr_accesses gt π π2 ->
    π1 = π2.
Proof.
  introv H1 H2. gen π2. induction H1; intros;
  inverts_head tr_accesses; repeat fequals*; 
  inverts_head access_field; subst; simpls;
  inverts_head Logic.or; repeat fequals*.
Qed. 

Theorem functional_tr_val : forall gt v v1 v2,
  tr_val gt v v1 ->
  tr_val gt v v2 ->
  v1 = v2.
Proof using.
  introv H1 H2. gen v2. induction H1; intros;
  inverts_head tr_val; fequals*; subst; simpls; tryfalse.
  { applys* functional_tr_accesses. }
  { applys* eq_of_extens. math. }
  { applys read_extens. 
    { inverts_head make_group_tr'. congruence. }
    { introv Hin. tests C: (i = fg).
      { inverts_head make_group_tr'.
        asserts_rewrite~ (s'0[fg0] = val_struct (typ_var Tg0) sg0).
        asserts_rewrite~ (s'[fg0] = val_struct (typ_var Tg0) sg).
        fequals. applys~ read_extens. introv Hk. 
        asserts_rewrite* (dom sg = dom sg0) in *. }
      { inverts_head make_group_tr'.
        asserts_rewrite~ (dom s' = dom s \- dom sg \u '{fg0}) in Hin.
        inverts Hin as Hin; tryfalse. inverts Hin as Hin Hnotin.
        asserts_rewrite* (dom sg = dom sg0) in *. } } }
  { applys read_extens.
    { congruence. }
    { introv Hin. 
      asserts_rewrite* (dom s' = dom s) in *. } }
Qed.

Lemma functional_tr_struct_op : forall gt op t1 tr1 tr2,
  is_struct_op op ->
  tr_struct_op gt (trm_app op (t1 :: nil)) tr1 ->
  tr_struct_op gt (trm_app op (t1 :: nil)) tr2 ->
  tr1 = tr2.
Proof.
  introv Hop Htr1 Htr2. induction Htr1; subst; 
  inverts H; inverts Htr2; inverts_head Logic.or;
  try inverts_head make_group_tr'; simpls;
  try inverts_head prim_struct_access;
  try inverts_head prim_struct_get;
  repeat fequals*.
Qed.

Theorem functional_tr_trm : forall gt t t1 t2,
  tr_trm gt t t1 ->
  tr_trm gt t t2 ->
  t1 = t2.
Proof.
  introv H1 H2. gen t2. induction H1; intros;
  try solve [ inverts H2 ; try subst ; repeat fequals* ].
  { inverts H2. fequals. applys* functional_tr_val. }
  { inverts_head tr_trm; subst.
    { forwards~: IHtr_trm t1'0. subst.
      applys* functional_tr_struct_op. }
    { false. } }
Qed.

Theorem functional_tr_stack_item : forall gt i i1 i2,
  tr_stack_item gt i i1 ->
  tr_stack_item gt i i2 ->
  i1 = i2.
Proof.
  introv Hi1 Hi2. 
  inverts Hi1 as H. inverts Hi2 as H'. 
  forwards*: functional_tr_val H H'.
Qed.

Theorem functional_tr_stack : forall gt S S1 S2,
  tr_stack gt S S1 ->
  tr_stack gt S S2 ->
  S1 = S2.
Proof.
  introv HS1 HS2. inverts HS1 as HS1. inverts HS2 as HS2. 
  gen S2. induction HS1; intros.
  { inverts~ HS2. }
  { inverts HS2 as HSy0 HS1'. fequals.
    { forwards*: functional_tr_stack_item H HSy0. }
    { applys~ IHHS1. } }
Qed.

Theorem functional_tr_state : forall gt m m1 m2,
  tr_state gt m m1 ->
  tr_state gt m m2 ->
  m1 = m2.
Proof.
  introv Hm1 Hm2. 
  inverts Hm1 as HD1 Htr1. inverts Hm2 as HD2 Htr2.
  applys read_extens.
  { unfolds state. congruence. }
  { introv Hi. rewrite HD1 in *. 
    forwards Hm1i: Htr1 Hi. forwards Hm2i: Htr2 Hi.
    forwards~: functional_tr_val Hm1i Hm2i. }
Qed.


(* ********************************************************************** *)
(* * Correctness of the transformation *)

Section TransformationsProofs.

(* ---------------------------------------------------------------------- *)
(** Path surgery *)

Lemma tr_accesses_app : forall gt π1 π2 π1' π2',
  tr_accesses gt π1 π1' ->
  tr_accesses gt π2 π2' ->
  tr_accesses gt (π1 ++ π2) (π1' ++ π2').
Proof.
  introv Ha1 Ha2. gen π2 π2'. induction Ha1; intros;
  rew_list in *; eauto. 
Qed.


(* ---------------------------------------------------------------------- *)
(** Regularity of the transformation with respect to values *)

Lemma not_is_val_tr : forall gt t1 t2,
  tr_trm gt t1 t2 ->
  ~ is_val t1 ->
  ~ is_val t2.
Proof.
  introv Htr Hv. induction Htr; introv HN;
  try solve [ subst ; inverts HN ]. forwards*: Hv.
  inverts_head tr_struct_op; inverts HN.
Qed.

Lemma not_is_error_tr : forall gt v1 v2,
  tr_val gt v1 v2 ->
  ~ is_error v1 ->
  ~ is_error v2.
Proof.
  introv Htr He. induction Htr; introv HN;
  try solve [ subst ; inverts HN ]. 
Qed.

Lemma not_is_uninitialized_tr : forall gt v v',
  tr_val gt v v' -> 
  ~ is_uninitialized v ->
  ~ is_uninitialized v'.
Proof.
  introv Htr Hu. induction Htr; introv HN;
  subst; inverts HN. forwards~: Hu.
Qed.


Lemma tr_accesses_inj : forall C gt π π1 π2,
  group_tr_ok gt C ->
  valid_accesses C π1 ->
  valid_accesses C π2 ->
  tr_accesses gt π1 π ->
  tr_accesses gt π2 π ->
    π1 = π2.
Proof.
  introv Hok Hva1 Hva2 Hπ1 Hπ2. gen C π2. induction Hπ1; intros.
  { inverts Hπ2. auto. }
  { inverts Hπ2; inverts Hva1; inverts Hva2.
    { fequals. applys* IHHπ1. }
    { fequals. } }
  { subst. inverts Hπ2; inverts Hva1; inverts Hva2.
    { fequals. applys* IHHπ1. }
    { simpls. inverts Hok as Hgt. inverts Hgt.
      inverts_head Logic.or; tryfalse. inverts H4.
      { inverts H14. fequals. }
      { inverts H8. inverts H12. false*. } } }
  { inverts Hπ2; inverts Hva1; inverts Hva2.
    { inverts Hok as Hgt. inverts Hgt. inverts H8.
      subst. simpls. inverts_head Logic.or; tryfalse.
      inverts H6. inverts Hπ1.
      { inverts H7. inverts H15. inverts H19. inverts H3; fequals. }
      { simpls. inverts H18.
        { inverts H7. inverts H15. false*. }
        { false*. } } }
    { fequals. applys* IHHπ1. } }
Qed.

(* TODO: Set lemmas *)
Lemma subset_union_eq : forall A (S1 S2:set A),
  S1 \c S2 ->
  S2 = (S2 \- S1) \u S1.
Proof.
  intros. rew_set in *. intuition.
  { rew_set. tests: (x \in S1).
    { right~. }
    { left~. } }
  { rew_set in *. inverts H0.
    { inverts~ H1. }
    { applys~ H. } }
Qed.

Lemma incl_eq : forall A (S1 S2 S:set A),
  S \c S1 ->
  S \c S2 ->
  (S1 \- S) = (S2 \- S) ->
  S1 = S2.
Proof.
  intros. 
  asserts_rewrite (S1 = (S1 \- S) \u S). 
  { applys~ subset_union_eq. }
  asserts_rewrite (S2 = (S2 \- S) \u S). 
  { applys~ subset_union_eq. }
  rew_set. intuition.
  { rewrite~ <- H1. }
  { rewrite~ H1. }
Qed.
(* ------------- *)

Lemma tr_val_inj : forall T1 T2 φ C gt v v1 v2,
  group_tr_ok gt C ->
  valid_phi C φ ->
  typing_val C φ v1 T1 ->
  typing_val C φ v2 T2 ->
  tr_val gt v1 v ->
  tr_val gt v2 v ->
  v1 = v2.
Proof.
  introv Hok Hφ HTv1 HTv2 Hv1 Hv2. gen C φ v2 T1 T2. induction Hv1; intros; 
  inverts Hv2; repeat fequals*; subst; simpls; tryfalse*.
  { inverts HTv1 as HRφ1. inverts HTv2 as HRφ2. 
    inverts HRφ1. inverts HRφ2. 
    unfolds valid_phi.
    forwards~: follow_typ_valid_accesses φ[l] T C π.
    forwards~: follow_typ_valid_accesses φ[l] T0 C π0.
    forwards~: Hφ l.
    applys* tr_accesses_inj. }
  { applys* eq_of_extens. congruence. 
    inverts HTv1. inverts HTv2. introv Hi.
    applys* H1. }
  { applys* read_extens. 
    { inverts_head make_group_tr'. rewrite <- H3 in *.
      applys~ incl_eq (dom s) (dom s0) (dom sg).
      rewrite H2 in H14. admit. (*True but more set lemmas needed.*) }
    { introv Hi. inverts HTv1. inverts HTv2.
      inverts_head make_group_tr'. 
      rewrite H8 in H19. inverts H19. tests: (i \indom sg0).
      { applys* H5.
        { applys~ H15. rewrite~ H10. }
        { applys~ H21. admit. admit. (*Set lemmas*) } }
      { applys* H7. 
        { applys~ H18. admit. (*More set lemmas needed*) } 
        { applys~ H15. rewrite~ H10. }
        { applys~ H21. admit. admit. (*Set lemmas*) } } } }
  { inverts HTv1. inverts HTv2. applys* read_extens. congruence.
    introv Hi. applys* H3. 
    { applys~ H10. admit. (*Some rewrites*) }
    { applys~ H11. rewrite~ H6. }
    { applys~ H14. admit. admit. (*Some rewrites*) } }
Qed.

Lemma tr_val_inj_cp : forall C φ T1 T2 gt v1 v2 v1' v2',
  group_tr_ok gt C ->
  valid_phi C φ ->
  typing_val C φ v1 T1 ->
  typing_val C φ v2 T2 ->
  tr_val gt v1 v1' ->
  tr_val gt v2 v2' ->
  v1 <> v2 ->
  v1' <> v2'.
Proof.
  introv Hok Hφ HTv1 HTv2 Hv1 Hv2 Hneq HN. subst. 
  forwards*: tr_val_inj Hok Hφ HTv1 HTv2 Hv1.
Qed.

Lemma not_tr_val_error : forall gt v1 v2,
  tr_val gt v1 v2 ->
  ~ is_error v2.
Proof.
  introv Hv He. unfolds is_error. destruct* v2. inverts Hv.
Qed.

Lemma tr_stack_add : forall gt z v S v' S',
  tr_stack gt S S' ->
  tr_val gt v v' ->
  tr_stack gt (Ctx.add z v S) (Ctx.add z v' S').
Proof.
  introv HS Hv. constructors~. inverts HS.
  unfolds Ctx.add. destruct* z.
  applys~ Forall2_cons. constructors~.
Qed.


(* ---------------------------------------------------------------------- *)
(** uninitialized is coherent with the transformation *)

Lemma tr_uninitialized_val' : forall gt v v' T C C',
  tr_typdefctx gt C C' ->
  tr_val gt v v' ->
  uninitialized C T v ->
  uninitialized C' T v'.
Proof using.
  introv HC Hv Hu. gen C' v'. induction Hu; intros;
  try solve [ inverts Hv ; constructors~ ].
  { (* val array *)
    inverts Hv as Hl Hai. constructors*. 
    2 : { rewrite <- Hl. eauto. }
    inverts_head typing_array.
    { constructors. }
    { inverts HC. constructors. admit. tests: (Tv=Tt).
      { inverts H4. 
        { rewrite H6 in H12. inverts H12. }
        { rewrite H6 in H. inverts H. } }
      { rewrite~ H9. inverts H4.
        { constructors. }
        { admit. } } } }
(*   { (* val struct *)
    tests: (T = group_tr_struct_name gt); 
    inverts Hv as; inverts HC as; 
    try solve [ intros ; simpls ; tryfalse ].
    { (* fields grouped *)
      introv HTtD HTgD HDC' HDCT HDCTt HfgD HCTtf HDC'Tt HC'Ttfg HC'Tgf.
      introv HDvfs Hmg Hfg0 HDs' Htrsg Htrs' Hs'fg0. 
      constructors~; unfolds typdefctx; unfolds typdef_struct.
      { rewrite HDC'. rew_set~. }
      { inverts Hmg. simpls. congruence. }
      { inverts Hmg. simpls. introv Hi. subst.
        tests: (f = fg0).
        { rewrite HC'Ttfg. rewrite Hs'fg0. constructors~. 
          { rewrite HDC' at 1. rew_set~. }
          { introv Hi'. rewrite~ HC'Tgf. forwards~: H3 f C' sg[f].
            { rew_set in *. applys~ HDCTt. }
            { constructors~. }
            { applys~ Htrsg. rewrite~ <- H5. } } }
        { forwards~: H3 f C' (s'[f]).
          { rewrite HDC'Tt in Hi; rew_set in *. 
            inverts Hi as (Hin&Hnin); eauto. }
          { constructors~. }
          { applys~ Htrs'; rewrite HDC'Tt in Hi; rew_set in Hi.
            { inverts Hi as (Hin&Hnin); tryfalse.
              rewrite~ <- H5. }
            { inverts Hi as (Hin&Hnin); tryfalse. 
              rewrite~ <- H1. } }
          { rewrite~ HCTtf; rewrite HDC'Tt in Hi; rew_set in Hi;
            inverts Hi as (Hin&Hnin); tryfalse; auto. } } } }
    { (* other struct *)
      introv HTtD HTgD HDC' HC'T HDC'Tg HfgD HC'Tt HDC'Tt HC'Ttfg HC'Tgf. 
      introv Hneq HDvfs Htrs'f.
      constructors~; unfolds typdefctx; unfolds typdef_struct.
      { rewrite HDC'. rew_set~. }
      { simpls. rewrite~ HC'T. congruence. }
      { introv Hi. rewrite~ HC'T. subst. forwards~: H3 f C' (s'[f]).
        { rewrite~ <- HC'T. }
        { constructors~. }
        { applys~ Htrs'f. rewrite <- H1. rewrite~ <- HC'T. } } } }
(* TODO: I had to add this because there were 'remaining goals on the shelf'. *)
Unshelve. typeclass. *) admit.
Qed.

(* This will be proved when the relation is translated to a function. 
   See TrGroupFun.v. *)
Lemma total_tr_val' : forall gt v,
  exists v', tr_val gt v v'.
Proof.
Admitted.

(* Lemma for the new case. *)
Lemma tr_uninitialized_val : forall gt v T C C',
  tr_typdefctx gt C C' ->
  uninitialized C T v ->
  exists v',
        tr_val gt v v'
    /\  uninitialized C' T v'.
Proof.
  introv HC Hu. forwards* (v'&Hv'): total_tr_val' gt v.
  exists v'. splits~. applys* tr_uninitialized_val'.
Qed.


(* ---------------------------------------------------------------------- *)
(** Correctness of access transformations *)

Lemma tr_read_accesses : forall gt v π v' π' w,
  tr_val gt v v' ->
  tr_accesses gt π π' ->
  read_accesses v π w ->
  (exists w',
      tr_val gt w w'
  /\  read_accesses v' π' w').
Proof.
  introv Hv Ha HR. gen gt v' π'. induction HR; intros.
  { (* nil *)
    inverts Ha. exists~ v'. }
  { (* array_access *)
    inverts Ha as Ha. inverts Hv as Hl Htr.
    forwards Htra: Htr H.
    forwards (w'&Hw'&Hπ'): IHHR Htra Ha.
    exists* w'. }
  { (* struct_access *)
    inverts Ha as; inverts Hv as;
    try solve [ intros ; false ].
    { (* one of the fields to group *)
      introv HD1 Hgt HD2 HD3 Hsg Hfs Hsfg Hπ Hin.
      rewrite Hgt in *. simpls.
      forwards Hsf: Hsg Hin.
      forwards (w'&Hw'&HR'): IHHR Hsf Hπ.
      exists w'. splits~.
      constructors. 
      { rewrite HD3. rew_set~. } 
      rewrite Hsfg. constructors~. }
    { (* struct transformed but another field *)
      introv HD1 HD2 HD3 Hsg Hfs Hsfg Hπ Hor.
      inverts Hor as Hf; simpl in Hf; tryfalse.
      forwards Hsf: Hfs Hf H.
      forwards (w'&Hw'&HR'): IHHR Hsf Hπ.
      exists w'. splits~. constructors~.
      rewrite HD3. rew_set~. }
    { (* another struct *)
      intros Hn HD Hfs Hπ Hor.
      forwards Hsf: Hfs H.
      forwards (w'&Hw'&HR'): IHHR Hsf Hπ.
      exists w'. splits~. constructors~.
      rewrite~ <- HD. } }
Qed.


Lemma tr_write_accesses : forall v1 w gt π v1' π' w' v2,
  tr_val gt v1 v1' ->
  tr_val gt w w' ->
  tr_accesses gt π π' ->
  write_accesses v1 π w v2 ->
  (exists v2',
        tr_val gt v2 v2'
    /\  write_accesses v1' π' w' v2').
Proof.
  introv Hv1 Hw Hπ HW. gen gt v1' w' π'. induction HW; intros.
  { (* nil *)
    inverts Hπ. exists~ w'. }
  { (* array_access *)
    inverts Hπ as Hπ. inverts Hv1 as Hl Htr.
    forwards Htra: Htr H.
    forwards (v2'&Hv2'&HW'): IHHW Htra Hw Hπ.
    exists (val_array a'[i:=v2']).
    splits; constructors*; rewrite H0.
    { repeat rewrite~ length_update. }
    { introv Hi0. rew_reads; rew_index* in *. } }
  { (* struct_access *)
    inverts Hπ as; inverts Hv1 as;
    try solve [ intros ; false ].
    { (* one of the fields to group *)
      introv HD1 Hgt HD2 HD3 Hsg Hfs Hsfg Hπ Hin.
      rewrite Hgt in *. simpls.
      forwards Hsf: Hsg Hin.
      forwards (v2'&Hv2'&HW'): IHHW Hsf Hw Hπ.
      remember (group_tr_struct_name gt) as T.
      exists (val_struct T s'[fg:=(val_struct Tsg sg[f:=v2'])]).
      substs. splits.
      { applys~ tr_val_struct_group (sg[f:=v2']);
        repeat rewrite dom_update_at_indom in *; eauto.
        { rewrite HD3. rew_set~. }
        { intros. rew_reads*. }
        { intros. rew_reads; intros; subst; eauto; contradiction. }
        { rew_reads~. } }
      { constructors~. rewrite HD3. rew_set~.
        rewrite Hsfg. constructors*. } }
    { (* struct transformed but another field *)
      introv HD1 HD2 HD3 Hsg Hfs Hsfg Hπ Hor.
      inverts Hor as Hf; simpl in Hf; tryfalse.
      forwards Hsf: Hfs Hf H.
      forwards (v2'&Hv2'&HW'): IHHW Hsf Hw Hπ.
      exists (val_struct T s'[f:=v2']). splits.
      { applys~ tr_val_struct_group; subst;
        repeat rewrite dom_update_at_indom in *; eauto.
        { rewrite HD3. rew_set~. }
        { intros. rew_reads; intros; eauto. 
          subst. forwards*: in_notin_neq. }
        { intros. rew_reads*. }
        { rew_reads~. intros. 
          subst. forwards*: in_notin_neq. } }
      { constructors~. rewrite HD3. rew_set~. auto. } }
    { (* another struct *)
      intros Hn HD Hfs Hπ Hor.
      forwards Hsf: Hfs H.
      forwards (v2'&Hv2'&HW'): IHHW Hsf Hw Hπ.
      exists (val_struct T s'[f:=v2']). splits.
      { constructors~; subst;
        repeat rewrite dom_update_at_indom; eauto.
        { rewrite~ <- HD. }
        { intros. rew_reads~. } }
      { constructors~. rewrite* <- HD. auto. } } }
Qed.


(* ---------------------------------------------------------------------- *)
(** Correctness of the transformation *)

(* TODO: Rewrites used quite often throughout: 
   - dom_update_at_indom
   - dom_update *)

(* TODO: Sometimes, in order to use map tactics or 
   lemmas, I have to do [unfold state] in order to get
   [map loc var] instead of [state] but then the
   rewrites stop working so I have to [fold state] again.
   How can we fix this? *)

Theorem red_tr: forall gt C C' t t' v S S' m1 m1' m2,
  tr_typdefctx gt C C' ->
  tr_trm gt t t' ->
  tr_stack gt S S' ->
  tr_state gt m1 m1' ->
  red C S m1 t m2 v ->
  ~ is_error v ->
  exists v' m2',
      tr_val gt v v'
  /\  tr_state gt m2 m2'
  /\  red C' S' m1' t' m2' v'.
Proof.
  introv HC Ht HS Hm1 HR He. gen C' t' S' m1'. induction HR; intros;
  try solve [ forwards*: He; unfolds* ].
  { (* var *)
    inverts Ht. forwards* (v'&H'&Hv'): stack_lookup_tr HS H.
    exists* v' m1'. }
  { (* val *)
    inverts Ht. exists* v' m1'. }
  { (* if *)
    inverts Ht as Hb HTrue HFalse.
    forwards* (v'&m2'&Hv'&Hm2'&HR3): IHHR1 Hb HS Hm1.
    introv HN. inverts HN.
    inverts* Hv'.
    forwards* (vr'&m3'&Hvr'&Hm3'&HR4): IHHR2 HS Hm2'. 
    2: exists* vr' m3'. case_if*. }
  { (* let *)
    inverts Ht as Ht1 Ht2.
    forwards* (v'&m2'&Hv'&Hm2'&HR3): IHHR1 Ht1 HS Hm1.
    forwards HS': tr_stack_add z HS Hv'.
    forwards* (vr'&m3'&Hvr'&Hm3'&HR4): IHHR2 Ht2 HS' Hm2'.
    forwards: not_tr_val_error Hv'.
    exists* vr' m3'. }
  { (* binop *)
    inverts Ht as Ht1 Ht2.
    inverts Ht1 as Ht1. inverts Ht2 as Ht2.
    inverts H.
    { exists __ m1'. splits~. inverts Ht1. inverts Ht2. 
      constructors. constructors~. }
    { exists __ m1'. splits~. inverts Ht1. inverts Ht2. 
      constructors. constructors~. }
    { exists __ m1'. splits~. constructors.
      forwards: functional_tr_val Ht1 Ht2. subst.
      constructors~; applys* not_is_error_tr. }
    { exists __ m1'. splits~. constructors.
      forwards: neq_tr gt H2 Ht1 Ht2.
      constructors~. (* TODO: try applys not_is_error_tr v1 v' *)
      { apply not_is_error_tr with (gt:=gt) (v1:=v1) (v2:=v'); auto. }
      { apply not_is_error_tr with (gt:=gt) (v1:=v2) (v2:=v'0); auto. } } }
  { (* get *)
    inverts Ht as _ Hp. subst.
    inverts Hm1 as HD Htrm.
    inverts H0 as Hi Ha.
    forwards Htrml: Htrm Hi.
    inverts Hp as Hp. inverts Hp as Hπ.
    forwards (w'&Hw'&Ha'): tr_read_accesses Htrml Hπ Ha.
    exists w' m1'. splits*.
    repeat constructors~. rewrite~ <- HD.
    applys* not_is_uninitialized_tr. }
  { (* set *)
    inverts Ht as Hp Ht. subst.
    inverts Hm1 as HD Htrm.
    inverts H2 as Hin HW. 
    forwards Htrml: Htrm Hin.
    inverts Hp as Hp. 
    inverts Hp as Hπ.
    inverts Ht as Hv.
    forwards (w'&Hw'&HW'): tr_write_accesses Htrml Hv Hπ HW.
    exists val_unit m1'[l:=w']. splits~.
    { constructors.
      { unfold state. repeat rewrite~ dom_update.
        fold state. rewrite~ HD. }
      { introv Hi'. rew_reads~. intros. applys Htrm.  
        applys~ indom_update_inv_neq Hi'. } }
    { constructors~. applys* not_tr_val_error.
      constructors*. rewrite~ <- HD. } }
  { (* new *) 
    inverts Ht. subst.
    inverts Hm1 as HD Htrm. 
    forwards* (v'&Hv'&Hu): tr_uninitialized_val.
    exists (val_abstract_ptr l0 nil) m1'[l0:=v']. splits~.    
    { constructors.
      { unfold state. repeat rewrite~ dom_update.
        fold state. rewrite~ HD. }
      { introv Hin. unfolds state. rew_reads; intros; eauto. } }
    { constructors~. rewrite~ <- HD. auto. } }  
  { (* new_array *)
    inverts Ht as _ Ht.
    inverts Ht as Hv.
    inverts Hv. 
    inverts Hm1 as HD Htrm. subst.
    forwards* (v''&Hv''&Hu): tr_uninitialized_val.
    exists (val_abstract_ptr l0 nil) m1'[l0:=v'']. splits~.
    { constructors.
      { unfold state. repeat rewrite~ dom_update.
        fold state. rewrite~ HD. }
      { introv Hin. unfolds state. rew_reads; intros; eauto. } }
    { applys~ red_new_array. rewrite~ <- HD. auto. } }
  { (* struct_access *)
    inverts Ht as; inverts Hm1 as HD Htrm.
    { (* accessing grouped field *)
      introv Ht Hf. subst.
      inverts Ht as Hv. inverts Hv as Ha.
      remember (group_tr_struct_name gt) as Ts.
      remember (group_tr_new_struct_name gt) as Tsg.
      remember (group_tr_new_struct_field gt) as fg.
      remember (access_field Ts fg) as a1.
      remember (access_field Tsg f) as a2.
      exists (val_abstract_ptr l (π'++(a1::a2::nil))) m1'.
      splits*.
      { constructors. applys* tr_accesses_app. subst*. }
      { subst. applys* red_args_1. applys* red_struct_access.
        fequals*. rew_list*. } }
    { (* accessing another field *) 
      introv Ht Hor. subst. inverts Ht as Hv. inverts Hv as Ha.
      exists (val_abstract_ptr l (π'++(access_field T f :: nil))) m1'.
      splits; constructors*. applys* tr_accesses_app. }
    { introv HN. forwards*: HN. } }
  { (* array_access *)
    inverts Ht as Ht Hti. subst.
    inverts Ht as Hv. inverts Hv as Ha.
    inverts Hti as Hv. inverts Hv.
    inverts Hm1 as HD Htrm.
    exists (val_abstract_ptr l (π'++(access_array i::nil))) m1'.
    splits; constructors*. applys* tr_accesses_app. }
  { (* struct_get *) 
    inverts Ht as.
    { (* accessing grouped field *)
      introv Ht Hf. subst.
      inverts Ht as Hv. inverts Hv as; 
      try solve [ intros ; contradiction ].
      introv HDsg Hgt Hfg HDs' Hsf Htrsf Hs'fg.
      exists sg[f] m1'.
      splits~.
      { applys~ Hsf. rewrite Hgt in *. simpls~. }
      { applys~ red_args_1.
        { applys~ red_struct_get. rewrite Hgt. simpls~.
          rewrite HDs'. rew_set~. }
        { rewrite Hgt in *. simpls. applys~ red_struct_get. 
          rewrite~ Hs'fg. } } }
    { (* accessing another field *) 
      introv Ht Hor. subst. inverts Ht as Hv. inverts Hor.  
      { inverts Hv as; try solve [ intros ; contradiction ].
        introv _ HDs Htrsf. exists s'[f] m1'. splits~. constructors~.
        rewrite~ <- HDs. }
      { inverts Hv as. 
        { introv HDsg Hfg HDs' Hsf Htrsf Hs'fg. 
          exists s'[f] m1'. splits~. constructors~.
          rewrite HDs'. rew_set~. }
        { introv Hneq HDs Htrsf. exists s'[f] m1'. splits~.
          constructors~. rewrite~ <- HDs. } } }
    { introv HN. forwards*: HN. } }
  { (* array_get *) 
    inverts Ht as Ht Hti. subst.
    inverts Ht as Hv.
    inverts Hti as Hvi.
    inverts Hv as Hl Hai.
    inverts Hvi. 
    exists a'[i] m1'. 
    splits~. constructors*. }
  { (* args_1 *)
    inverts Ht; forwards* (v'&m2'&Hv'&Hm2'&HR'): IHHR1;
    forwards*: not_is_error_args_1 HR2 He;
    forwards* (v''&m3'&Hv''&Hm3'&HR''): IHHR2;
    exists v'' m3'; splits*;
    try solve [ applys* red_args_1; applys* not_is_val_tr ].
    { (* Case struct access *)
      inverts HR''.
      { tryfalse. }
      { applys* red_args_1. applys* red_args_1.
        apply not_is_val_tr with (gt:=gt) (t1:=t1); auto. }
      { forwards HN: not_is_error_tr He Hv''. forwards*: HN. unfolds*. }
      { forwards HN: not_is_error_tr He Hv''. forwards*: HN. unfolds*. } }
    { (* Case struct get *) 
      inverts HR''. 
      { tryfalse. } 
      { applys* red_args_1. applys* red_args_1.
        apply not_is_val_tr with (gt:=gt) (t1:=t1); auto. }
      { forwards HN: not_is_error_tr He Hv''. forwards*: HN. unfolds*. } } }
  { (* args_2 *)
    inverts Ht as Ht1 Ht2.
    forwards* (v'&m2'&Hv'&Hm2'&HR'): IHHR1.
    forwards*: not_is_error_args_2 HR2 He.
    forwards* (v''&m3'&Hv''&Hm3'&HR''): IHHR2.
    exists v'' m3'. splits*.
    inverts Ht1. applys* red_args_2.
    applys* not_is_val_tr. }
Qed.

End TransformationsProofs.
