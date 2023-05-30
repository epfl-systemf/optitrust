open Ast
open Stats
open Tools
open PPrint

(* [ml_file_excerpts]: maps line numbers to the corresponding sections in-between [!!] marks in
   the source file. Line numbers are counted from 1 in that map. *)
module Int_map = Map.Make(Int)
let ml_file_excerpts = ref Int_map.empty

(* [compute_ml_file_excerpts lines]: is a function for grouping lines according to the [!!] symbols. *)
let compute_ml_file_excerpts (lines : string list) : string Int_map.t =
  let r = ref Int_map.empty in
  let start = ref 0 in
  let acc = Buffer.create 3000 in
  let push () =
    r := Int_map.add (!start+1) (Buffer.contents acc) !r;
    Buffer.clear acc; in
  let regexp_let = Str.regexp "^[ ]*let" in
  let starts_with_let (str : string) : bool =
    Str.string_match regexp_let str 0 in
  (* match a line that starts with '!!' or 'bigstep' *)
  let regexp_step = Str.regexp "^[ ]*\\(!!\\|bigstep\\)" in
  let starts_with_step (str : string) : bool =
    Str.string_match regexp_step str 0 in
  let process_line (iline : int) (line : string) : unit =
    if starts_with_step line then begin
      push();
      start := iline;
    end;
    if not (starts_with_let line) then begin
      Buffer.add_string acc line;
      Buffer.add_string acc "\n";
    end;
    in
  List.iteri process_line lines;
  push();
  !r


(******************************************************************************)
(*                             Logging management                             *)
(******************************************************************************)


(* [now()] returns the current time *)
let now () : float =
   Unix.gettimeofday()

(* [timing_log_handle]: is a handle on the channel for writing timing reports. *)
let timing_log_handle = ref None

(* [stats_log_handle]: is a handle on the channel for writing stats reports. *)
let stats_log_handle = ref None

(* [logs]: is a reference on the list of open log channels. *)
let logs : (out_channel list) ref = ref []

(* [close_logs]: closes all open log channels. *)
let close_logs () : unit =
  List.iter (fun log -> close_out log) !logs;
  logs := []

(* [init_logs]: initializes the log files. It closes any existing logs.
   Returns the one log created. *)
let init_logs prefix =
  close_logs();
  let clog = open_out (prefix ^ ".log") in
  let timing_log = open_out ("timing.log") in
  timing_log_handle := Some timing_log;
  let stats_log = open_out ("stats.log") in
  stats_log_handle := Some stats_log;
  logs := timing_log :: stats_log :: clog :: [];
  clog

(* [write_log clog msg]: writes the string [msg] to the channel [clog]. *)
let write_log (clog : out_channel) (msg : string) : unit =
  output_string clog msg;
  flush clog

(* [trm_to_log clog styp t]: writes in the channel [clog] the term [t],
   and its typing information described by the string [styp]. *)
let trm_to_log (clog : out_channel) (exp_type : string) (t : trm) : unit =
  let sloc =
    match t.loc with
    | None -> ""
    | Some {loc_file = _; loc_start = {pos_line = start_row; pos_col = start_column}; loc_end = {pos_line = end_row; pos_col = end_column}} ->
       Printf.sprintf "at start_location %d  %d end location %d %d" start_row start_column end_row end_column
    in
  let msg = Printf.sprintf " -expression\n%s\n %s is a %s\n" (AstC_to_c.ast_to_string t) sloc exp_type in
 write_log clog msg

(******************************************************************************)
(*                             File input                                     *)
(******************************************************************************)

(* A parser should read a filename and return:
   - A header to copy in the produced file (typically a list of '#include' for C)
   - The OptiTrust AST of the rest of the file *)
(* TODO: encode header information in the AST *)
type parser = string -> string * trm

(* [parse ~parser filename]:
   call the parser on the given file while recording statistics *)
let parse ~(parser: parser) (filename : string) : string * trm =
  print_info None "Parsing %s...\n" filename;
  let parsed_file = stats ~name:"tr_ast" (fun () -> parser filename) in
  print_info None "Parsing Done.\n";
  parsed_file

(******************************************************************************)
(*                             Trace management                               *)
(******************************************************************************)

(* [context]: contains general information about:
   - the source code that was loaded initially using [set_init_file],
   - the prefix of the filenames in which to output the final result using [dump]
   - the log file to report on the transformation performed. *)
type context =
  { parser : parser;
    (* DEPRECATED?
       directory : string; *)
    mutable prefix : string; (* TODO: needs mutable? *)
    extension : string;
    header : string;
    clog : out_channel; }

(* [contex_dummy]: used for [trace_dummy]. *)
let context_dummy : context =
  { parser = (fun _ -> failwith "context_dummy has no parser");
    (* directory = ""; *)
    prefix = "";
    extension = "";
    header = "";
    clog = stdout; }

(* DEPRECATED [stepdescr]: description of a script step. *)
type stepdescr = {
  mutable isbigstep : string option; (* if the step is the beginning of a big step,
                                        then the value is [Some descr] *)
  mutable script : string; (* excerpt from the transformation script, or "" *)
  mutable exectime : int; } (* number of milliseconds, -1 if unknown *)



(* [step_kind] : classifies the kind of steps *)
type step_kind = Step_root | Step_big | Step_small | Step_transfo | Step_target_resolve | Step_parsing

(* [step_kind_to_string] converts a step-kind into a string *)
let step_kind_to_string (k:step_kind) : string =
  match k with
  | Step_root -> "Root"
  | Step_big -> "Big"
  | Step_small -> "Small"
  | Step_transfo -> "Transfo"
  | Step_target_resolve -> "Target"
  | Step_parsing -> "Parsing"

(* [step_infos] *)
type step_infos = {
  mutable step_script : string;
  mutable step_script_line : int option;
  mutable step_time_start : float; (* seconds since start *)
  mutable step_duration : float; (* seconds *)
  mutable step_name : string;
  mutable step_args : (string * string) list;
  mutable step_justif : string list; }

(* [step_tree]: history type used for storing all the trace information about all steps, recursively *)
type step_tree = {
  mutable step_kind : step_kind;
  mutable step_ast_before : trm;
  mutable step_ast_after : trm;
  mutable step_sub : step_tree list;
  (* substeps in reverse order during construction (between open and close) *)
  mutable step_infos : step_infos; }
  (* TODO: FOR PARSING STEPS int_of_float(stats_parse.stats_time); } *)


(* A [step_stack] is a stack that contains the currently opened steps,
   with the innermost at the top. The bottom element of the stack is
   always the one that describes the execution of the full transformation
   script. *)
type step_stack = step_tree list


(* [trace]: a record made of a context, a current AST, and a list of ASTs that were
   saved as "interesting intermediate steps", via the [Trace.save] function.
   Any call to the [step] function adds a copy of [cur_ast] into [history]. *)
type trace = {
  mutable context : context;
  mutable cur_ast : trm;
  mutable step_stack : step_stack; } (* stack of open steps *)

(* [trm_dummy]: dummy trm. *)
let trm_dummy : trm =
  trm_val (Val_lit Lit_unit)

(* [trace_dummy]: an trace made of dummy context and dummy trm,
   whose purpose is to enforce that [Trace.init] is called before any transformation *)
let trace_dummy : trace =
  { context = context_dummy;
    cur_ast = trm_dummy; (* dummy *)
    step_stack = []; (* dummy *)
    }

(* [the_trace]: the trace produced by the current script. *)
let the_trace : trace =
  trace_dummy

(* [is_trace_dummy()]: returns whether the trace was never initialized. *)
let is_trace_dummy () : bool =
  the_trace.context == context_dummy

(* [get_decorated_history]: gets history from trace with a few meta information *)
(* TODO: remove this function? *)
let get_decorated_history ?(prefix : string = "") () : string * context * step_tree =
  let ctx = the_trace.context in
  let prefix = (* LATER: figure out why we needed a custom prefix here *)
    if prefix = "" then ctx.prefix else prefix in
  let tree =
    match the_trace.step_stack with
    | [] -> failwith "step stack must never be empty"
    | [tree] -> tree
    | _ -> failwith "step stack contains more than one step; this should not be the case when a transformation script has completed"
    in
  (prefix, ctx, tree)

let dummy_duration : float = 0.

(* [get_cur_step ()] returns the current step --there should always be one. *)
let get_cur_step () : step_tree =
  match the_trace.step_stack with
  | [] -> failwith "Trace.init has not been called"
  | step::_ -> step

(* [open_root_step] is called only by [Trace.init], for initializing
   the bottom element of the [step_stack].
   Assumes fields of [the_trace] have already been initialized. *)
let open_root_step ?(source : string = "<unnamed-file>") () : unit =
  assert(the_trace.step_stack = []);
  let step_root_infos = {
    step_script = "Contents of " ^ source;
    step_script_line = None;
    step_time_start = now();
    step_duration = dummy_duration;
    step_name = "Full script";
    step_args = [("extension", the_trace.context.extension) ];
    step_justif = [] }
   in
  let step_root = {
    step_kind = Step_root;
    step_ast_before = the_trace.cur_ast;
    step_ast_after = trm_dummy;
    step_sub = [];
    step_infos = step_root_infos; }
    in
  the_trace.step_stack <- [step_root]

(* [finalize_step] is called by [close_root_step] and [close_step] *)
let finalize_step (step : step_tree) : unit =
  let infos = step.step_infos in
  infos.step_duration <- now() -. infos.step_time_start;
  step.step_ast_after <- the_trace.cur_ast;
  step.step_sub <- List.rev step.step_sub

(* [get_root_step()] returns the root step, after close_root_step has been called *)
let get_root_step () : step_tree =
  match the_trace.step_stack with
  | [step] ->
      if step.step_ast_after == trm_dummy
        then failwith "get_root_step: close_root_step has not been called";
      step
  | _ -> failwith "close_root_step: broken invariant, stack must have size one"

(* [get_excerpt line]: returns the piece of transformation script that starts on the given line. Currently returns the ""
    in case [compute_ml_file_excerpts] was never called. LATER: make it fail in that case. *)
let get_excerpt (line : int) : string =
  if line = - 1 then "" else (*failwith "get_excerpt: requires a valid line number";*)
  if !ml_file_excerpts = Int_map.empty then "" else begin (* should "" be failure? *)
  match Int_map.find_opt line !ml_file_excerpts with
    | Some txt -> txt
    | None -> (*LATER: failwith? *) Printf.sprintf "<unable to retrieve line %d from script>" line
  end

(* [open_step] is called at the start of every big-step, or small-step,
   or combi, or basic transformation. *)
let open_step ?(line : int option) ?(step_script:string="") ~(kind:step_kind) ~(name:string) () : unit =
  let infos = {
    step_script;
    step_script_line = line;
    step_time_start = now();
    step_duration = dummy_duration;
    step_name = name;
    step_args = [];
    step_justif = [] }
   in
  let step = {
    step_kind = kind;
    step_ast_before = the_trace.cur_ast;
    step_ast_after = trm_dummy;
    step_sub = [];
    step_infos = infos; }
    in
  the_trace.step_stack <- step :: the_trace.step_stack

(* [step_justif txt] is called by a transformation after open_step in order
   to store an textual explaination of why it is correct. *)
let step_justif (justif:string) : unit =
  let step = get_cur_step() in
  let infos = step.step_infos in
  infos.step_justif <- justif::infos.step_justif

(* [step_justif_always_correct()] is a specialized version of [step_justif]
   for transformation that are always correct. *)
let step_justif_always_correct () : unit =
  step_justif "Transformation always correct"

(* [step_arg] is called by a transformation after open_step in order
   to store the string representations of one argument. *)
let step_arg ~(name:string) ~(value:string) : unit =
  let step = get_cur_step() in
  let infos = step.step_infos in
  infos.step_args <- (name,value)::infos.step_args

(* [close_step] is called at the end of every big-step, or small-step,
   or combi, or basic transformation. *)
let close_step () : unit =
  match the_trace.step_stack with
  | [] -> failwith "close_step: the_trace should not be empty"
  | [root_step] -> failwith "close_step: on the root, should call close_root_step"
  | step :: ((parent_step :: _) as stack_tail)  ->
      finalize_step step;
      parent_step.step_sub <- step :: parent_step.step_sub;
      the_trace.step_stack <- stack_tail

(* [close_step_kind_if_needed k] is used by
   [close_smallstep_if_needed] and [close_bigstep_if_needed] *)
let close_step_kind_if_needed (k:step_kind) : unit =
  let step = get_cur_step() in
  if step.step_kind = k then close_step()

(* [close_smallstep_if_needed()] closes a current big-step.
   Because big-steps are not syntactically scoped in the user script,
   we need such an implicit close operation to be called on either
   the opening of a new big-step, or on closing of the root step. *)
let close_smallstep_if_needed () : unit =
  close_step_kind_if_needed Step_small

(* [close_bigstep_if_needed()] closes a current big-step.
   Because big-steps are not syntactically scoped in the user script,
   we need such an implicit close operation to be called on either
   the opening of a new big-step, or on closing of the root step. *)
let close_bigstep_if_needed () : unit =
  close_smallstep_if_needed();
  close_step_kind_if_needed Step_big

(* [close_root_step] is called only by [Run. TODO???] at the end of the script,
   or at the end of the small-step targeted by the user.
   It leaves the root step at the bottom of the stack *)
let close_root_step () : unit =
  close_bigstep_if_needed();
  let step = match the_trace.step_stack with
    | [step] -> step
    | _ -> failwith "close_root_step: broken invariant, stack must have size one" in
  finalize_step step

(* [step] is a function wrapping the body of a transformation *)
let step ?(line : int = -1) ~(kind:step_kind) ~(name:string) (body : unit -> 'a) : 'a =
  open_step ~line ~kind ~name ();
  let r = body() in
  close_step();
  r

(* [parsing_step f] accounts for a parsing operation *)
let parsing_step (f : unit -> unit) : unit =
  step ~kind:Step_parsing ~name:"" f

(* [open_target_resolve_step] *)
let open_target_resolve_step () : unit =
  open_step ~kind:Step_target_resolve ~name:"" ()

(* [close_target_resolve_step] has a special handling because it saves a diff
   between an AST and an AST decorated with marks for targeted paths,
   even though the [cur_ast] is not updated with the marsk *)
let close_target_resolve_step (ps:Path.path list) (t:trm) : unit =
  if !Flags.dump_trace then begin
    let marked_ast, _marks = Path.add_marks_at_paths ps t in
    let cur_ast = the_trace.cur_ast in
    the_trace.cur_ast <- marked_ast;
    close_step();
    the_trace.cur_ast <- cur_ast
  end else
    close_step()

(* [invalidate()]: restores the global state (object [trace]) in its uninitialized state,
   like at the start of the program.  *)
let invalidate () : unit =
  close_logs();
  the_trace.context <- trace_dummy.context;
  the_trace.cur_ast <- trace_dummy.cur_ast;
  the_trace.step_stack <- trace_dummy.step_stack

(* [get_initial_ast ~parser ser_mode ser_file filename]: gets the initial ast before applying any trasformations
     [parser] - choose which parser to use for parsing the source code
     [ser_mode] - serialization mode
     [ser_file] - if serialization is used for the initial ast, the filename of the serialized version
                  of the source code is needed
     [filename] - filename of the source code  *)
let get_initial_ast ~(parser : parser) (ser_mode : Flags.serialization_mode) (ser_file : string)
  (filename : string) : (string * trm) =
  (* LATER if ser_mode = Serialized_Make then let _ = Sys.command ("make " ^ ser_file) in (); *)
  let ser_file_exists = Sys.file_exists ser_file in
  let ser_file_more_recent = if (not ser_file_exists) then false else Xfile.is_newer_than ser_file filename in
  let auto_use_ser = (ser_mode = Serialized_Auto && ser_file_more_recent) in
  if (ser_mode = Serialized_Use (* || ser_mode = Serialized_Make *) || auto_use_ser) then (
    if not ser_file_exists
      then fail None "Trace.get_initial_ast: please generate a serialized file first";
    if not ser_file_more_recent
      then fail None (Printf.sprintf "Trace.get_initial_ast: serialized file is out of date with respect to %s\n" filename);
    let ast = Xfile.unserialize_from ser_file in
    if auto_use_ser
      then Printf.printf "Loaded ast from %s.\n" ser_file;
    ast
  )
  else
    parse ~parser filename

(* [init f]: initializes the trace with the contents of the file [f].
   This operation should be the first in a transformation script.
   The history is initialized with the initial AST.
   [~prefix:"foo"] allows to use a custom prefix for all output files,
   instead of the basename of [f]. *)
(* LATER for mli: val set_init_source : string -> unit *)
let init ?(prefix : string = "") ~(parser: parser) (filename : string) : unit =
  invalidate ();
  let basename = Filename.basename filename in
  let extension = Filename.extension basename in
  let default_prefix = Filename.remove_extension filename in
  let ml_file_name =
    if Tools.pattern_matches "_inlined" default_prefix
      then List.nth (Str.split (Str.regexp "_inlined") default_prefix) 0
      else default_prefix in
  if !Flags.analyse_stats || !Flags.dump_trace then begin
    let src_file = (ml_file_name ^ ".ml") in
    if Sys.file_exists src_file then begin
      let lines = Xfile.get_lines src_file in
      (* printf "%s\n" (Tools.list_to_string ~sep:"\n" lines); *)
      ml_file_excerpts := compute_ml_file_excerpts lines;
    end;
  end;
  let mode = !Flags.serialization_mode in
  start_stats := get_cur_stats ();
  last_stats := !start_stats;

  let prefix = if prefix = "" then default_prefix else prefix in
  let clog = init_logs prefix in
  let ser_file = basename ^ ".ser" in

  let (header, cur_ast), stats_parse = Stats.measure_stats (fun () -> get_initial_ast ~parser mode ser_file filename) in

  let context = { parser; extension; prefix; header; clog } in
  the_trace.context <- context;
  the_trace.cur_ast <- cur_ast;
  the_trace.step_stack <- [];
  open_root_step ~source:ml_file_name ();

  if mode = Serialized_Build || mode = Serialized_Auto
    then Xfile.serialize_to ser_file (header, cur_ast);
  if mode = Serialized_Build
    then exit 0;
  print_info None "Starting script execution...\n"

(* [finalize()]: should be called at the end of the script to close the root step *)
let finalize () : unit =
  close_root_step()

(* [alternative f]: executes the script [f] in the original state that
   was available just after the call to [init].
   After the call, all the actions performed are discarded.

  Current usage:
     !! Trace.alternative (fun () ->
        !! Loop.fusion_on_block [cLabel "tofusion"];
        !!());

   LATER: figure out if it is possible to avoid "!!" in front and tail of [Trace.restart].
   LATER: figure out if this implementation could be extended in the presence of [switch]. *)
let alternative (f : unit->unit) : unit =
  failwith "unimplemented"
  (* TODO: fix this
  let trace = the_trace in
  if trace.history = [] || trace.stepdescrs = []
    then fail None "Trace.alternative: the history is empty";
  let _,init_ast = Xlist.unlast trace.history in
  let _,init_stepdescr = Xlist.unlast trace.stepdescrs in
  let cur_ast = trace.cur_ast in
  let history = trace.history in
  let stepdescrs = trace.stepdescrs in
  the_trace.cur_ast <- init_ast;
  the_trace.history <- [init_ast];
  the_trace.stepdescrs <- [init_stepdescr];
  f();
  the_trace.cur_ast <- cur_ast;
  the_trace.history <- history;
  the_trace.stepdescrs <- stepdescrs
  *)
  (* TODO: beautify? *)

(* [switch cases]: allows to introduce a branching point in a script.
   The [cases] argument gives a list of possible continuations (branches).
   Each of the branches can terminate with a [dump] operation to produce
   its output in a specific file. Alternatively, there could be further
   tranformations after the [switch] construct---a diamond construct.
   In such case, the instructions that follow the [switch] are applied
   in parallel on each of the traces, where one trace corresponds to one
   possible path in the script (via the branches).
   The optional argument [only_branch] can be use to temporary disable
   all branches but one. This is currently needed for the interactive mode
   to work. Branches are numbered from 1 (not from zero). *)
(* LATER for mli: switch : ?only_branch:int -> (unit -> unit) list -> unit *)
(* DEPRECATED
let switch ?(only_branch : int = 0) (cases : (unit -> unit) list) : unit =
  (* Close logs: new logs will be opened in every branch. *)
  close_logs ();
  let list_of_traces =
    Xlist.fold_lefti
      (fun i tr f ->
        let branch_id = i + 1 in
        if only_branch = 0 || branch_id = only_branch then
          begin
            let old_traces = !traces in
            let new_traces =
              List.fold_right
                (fun trace acc_traces ->
                  let context = trace.context in
                  (* create an extended prefix for this branch, unless there is a single branch *)
                  let prefix =
                    if List.length cases <= 1 || only_branch <> 0
                      then context.prefix
                      else context.prefix ^ "_" ^ (string_of_int branch_id)
                    in
                  (* create and register new log channel *)
                  let clog = open_out (context.directory ^ prefix ^ ".log") in
                  logs := clog :: !logs;
                  (* execute each branch in a single context *)
                  let branch_trace = { trace with context = { context with prefix; clog } } in
                  traces := [branch_trace];
                  f ();
                  (* store the traces produced by this branch *)
                  (!traces) :: acc_traces;
                )
                old_traces
                []
            in
            traces := old_traces;
            (List.flatten new_traces) :: tr
          end
        else tr
      )
      []
      cases
  in
  traces := List.flatten (List.rev list_of_traces)
 *)

(* [apply f]: applies the transformation [f] to the current AST,
   and updates the current ast with the result of that transformation.
   If there are several active trace (e.g., after a [switch]),
   then [f] is applied to each of the traces. During the execution of [f]
   on a given trace, the set of traces is replaced with a singleton set
   made of only that trace; this allows for safe re-entrant calls
   (i.e., the function [f] itself may call [Trace.apply]. *)
let apply (f : trm -> trm) : unit =
  if is_trace_dummy()
    then fail None "Trace.init must be called prior to any transformation.";
  the_trace.cur_ast <- f the_trace.cur_ast

(* [call f]: is similar to [apply] except that it applies to a function [f]
   with unit return type: [f] is meant to update the [cur_ast] by itself
   through calls to [apply].
   If there are several active trace (e.g., after a [switch]),
   then [f] is applied to each of the traces. During the execution of [f]
   on a given trace, the set of traces is replaced with a singleton set
   made of only that trace; this allows for safe re-entrant calls
   (i.e., the function [f] itself may call [Trace.apply]. *)
   (* TODO: see whether it's not simpler to use Trace.get_ast() ; DEPRECATED? *)
let call (f : trm -> unit) : unit =
  if is_trace_dummy()
    then fail None "Trace.init must be called prior to any transformation.";
  f the_trace.cur_ast




(******************************************************************************)
(*                                   Output                                   *)
(******************************************************************************)

(* [cleanup_cpp_file_using_clang_format filename]: makes a system call to
   reformat a CPP file using the clang format tool.
   LATER: find a way to remove extra parentheses in ast_to_doc, by using
   priorities to determine when parentheses are required. *)
let cleanup_cpp_file_using_clang_format ?(uncomment_pragma : bool = false) (filename : string) : unit =
  stats ~name:(Printf.sprintf "cleanup_cpp_file_using_clang_format(%s)" filename) (fun () ->
    ignore (Sys.command ("clang-format -style=\"Google\" -i " ^ filename));
    if (* temporary *) false && uncomment_pragma
      then ignore (Sys.command ("sed -i 's@//#pragma@#pragma@' " ^ filename))
  )


(* [get_header ()]: get the header of the current file (e.g. include directives) *)
let get_header () : string =
  the_trace.context.header

(* [output_prog ctx prefix ast]: writes the program described by the term [ast]
   in several files:
   - one describing the raw AST ("prefix.ast")
   - one describing the internal AST ("prefix_enc.cpp")
   - one describing the CPP code ("prefix.cpp").
   The CPP code is automatically formatted using clang-format. *)
let output_prog ?(beautify:bool=true) ?(ast_and_enc:bool=true) (ctx : context) (prefix : string) (ast : trm) : unit =
  let use_clang_format = beautify && !Flags.use_clang_format in
  let file_prog = prefix ^ ctx.extension in
  let out_prog = open_out file_prog in
  begin try
    (* print C++ code with decoding *)
    (*   DEPRECATED
    Printf.printf "===> %s \n" (ctx.includes); print_newline();*)
    (* LATER: try to find a way to put the includes in the AST so we can do simply ast_to_file *)
    output_string out_prog ctx.header;
    let beautify_mindex = beautify && !Flags.pretty_matrix_notation in
    if !Flags.bypass_cfeatures
      then AstC_to_c.ast_to_outchannel ~optitrust_syntax:true out_prog ast
      else AstC_to_c.ast_to_outchannel ~beautify_mindex ~comment_pragma:use_clang_format out_prog (Ast_fromto_AstC.cfeatures_intro ast);
    output_string out_prog "\n";
    close_out out_prog;
  with | Failure s ->
    close_out out_prog;
    failwith s
  end;
  (* beautify the C++ code --comment out for debug *)
  if use_clang_format
    then cleanup_cpp_file_using_clang_format ~uncomment_pragma:use_clang_format file_prog;
  (* ast and enc *)
  if ast_and_enc && !Flags.dump_ast_details then begin
    let file_ast = prefix ^ ".ast" in
    let file_enc = prefix ^ "_enc" ^ ctx.extension in
    let out_ast = open_out file_ast in
    let out_enc = open_out file_enc in
    begin try
      (* print the raw ast *)
      begin
        Ast_to_text.print_ast out_ast ast;
        output_string out_ast "\n";
        close_out out_ast;
      end;
      (* print the non-decoded ast *)
      output_string out_enc ctx.header;
      AstC_to_c.ast_to_outchannel ~optitrust_syntax:true out_enc ast;
      output_string out_enc "\n";
      close_out out_enc;
      if use_clang_format
        then cleanup_cpp_file_using_clang_format file_enc;
    with | Failure s ->
      close_out out_ast;
      close_out out_enc;
      failwith s
    end
  end

(******************************************************************************)
(*                                   Reparse                                  *)
(******************************************************************************)

(* [reparse_trm ctx ast]: prints [ast] in a temporary file and reparses it using Clang. *)
let reparse_trm ?(info : string = "") ?(parser: parser option) (ctx : context) (ast : trm) : trm =
  if !Flags.debug_reparse then begin
    let info = if info <> "" then info else "of a term during the step starting at" in
    Printf.printf "Reparse: %s.\n" info;
    flush stdout
  end;
  let in_prefix = (Filename.dirname ctx.prefix) ^ "/tmp_" ^ (Filename.basename ctx.prefix) in
  output_prog ~beautify:false ctx in_prefix ast;

  let parser =
    match parser with
    | Some p -> p
    | None -> ctx.parser
  in

  let (_, t) = parse ~parser (in_prefix ^ ctx.extension) in
  (*let _ = Sys.command ("rm " ^ in_prefix ^ "*") in*)
  t

(* [reparse ()]: function takes the current AST, prints it to a file, and parses it
   as if it was a fresh input. Doing so ensures in particular that all the type
   information is properly set up. WARNING: reparsing discards all the marks in the AST. *)
let reparse ?(info : string = "") ?(parser: parser option) () : unit =
  parsing_step (fun () ->
    let info = if info <> "" then info else "the code during the step starting at" in
    the_trace.cur_ast <- reparse_trm ~info ?parser the_trace.context the_trace.cur_ast
  )

(* Work-around for a name clash *)
let reparse_alias = reparse


(******************************************************************************)
(*                                   More dump                                *)
(******************************************************************************)


(* [dump_steps]: writes into files called [`prefix`_$i_out.cpp] the contents of each of the big steps,
    where [$i] denotes the index of a big step. *)
let dump_steps ?(onlybig : bool = false) ?(prefix : string = "") (foldername : string) : unit =
  ()
  (* TODO FIXME
  ignore (Sys.command ("mkdir -p " ^ foldername));
  let (prefix, ctx, hist_and_descr) = get_decorated_history ~prefix () in

  (* TODO: modify this code *)
  let n = List.length hist_and_descr in
  let id = ref 0 in
  List.iteri (fun i (ast,stepdescr) ->
    let isstartofbigstep =
      match stepdescr.isbigstep with
      | None -> false
      | Some _descr -> true
      in
    let should_dump = if onlybig then ((isstartofbigstep) || (i = n-1)) else true in
    if should_dump then begin
      let prefixi = Printf.sprintf "%s/%s_%s%d_out" foldername prefix (if !id < 10 then "0" else "") !id in
      output_prog ctx prefixi ast;
      incr id;
    end;
  ) hist_and_descr
  *)


(* [dump_trace_to_js]: writes into a file called [`prefix`.js] the
   contents of each of the steps record by the script, both for
   small steps and big steps, including the diffs and the excerpts
   of proof scripts associated with each step.

   The argument [history_and_isbigstep] takes the history
   with oldest entry first (unliked the history record field).
   It does not include the parsing step. --LATER: include it.

   The JS file is
   structured as follows (up to the order of the definitions):

   var codes = []; // if the script has 3 '!!', the array will have 4 entries (one for the state of the code before each '!!' and one for the final result)
   codes[i] = window.atob("...");

   var smallsteps = []; // smallsteps.length = codes.length - 1
   smallsteps[i] = { diff: window.atob("...");
                     script: window.atob("...");
                     exectime: ... } // for future use
   var bigsteps = []; // bigsteps.length <= smallsteps.length
   bigstep.push ({ diff: window.atob("...");
                  start: idStart;
                  stop: idStop;
                  descr : window.atob("...") });
     // invariant: bigstep[j].stop = bigstep[j+1].start
   *)
   (* TODO NEW

    step[i] = {  ... infos ; sub : [ i1; i2 ; ... ] }

    unique id for each step node in the tree

   *)
     (*
let dump_trace_to_js (history : history) : unit =
  ()

  let (prefix, ctx, hist_and_descr) = history in
  let file_js = prefix ^ "_trace.js" in
  let out_js = open_out file_js in
  let out = output_string out_js in
  let sprintf = Printf.sprintf in
  let cmd s =
    (* FOR DEBUG Printf.printf "execute: %s\n" s; flush stdout; *)
    ignore (Sys.command s) in
  let compute_command_base64 (s : string) : string =
    cmd (sprintf "%s | base64 -w 0 > tmp.base64" s);
    Xfile.get_contents ~newline_at_end:false "tmp.base64"
    in
  let compute_diff () : string =
    compute_command_base64 "git diff --ignore-all-space --no-index -U10 tmp_before.cpp tmp_after.cpp" in
  let lastbigstepstart = ref (-1) in
  let nextbigstep_descr = ref "<undefined bigstep descr>" in
  (* LATER: catch failures *)
  (* LATER: support other languages than C/C++ *)
  out "var codes = [];\nvar smallsteps = [];\nvar bigsteps = [];\n";
  let n = List.length hist_and_descr in
  List.iteri (fun i (ast,stepdescr) ->
    (* obtain source code *)
    output_prog ctx "tmp_after" ast;
    let src = compute_command_base64 "cat tmp_after.cpp" in
    out (sprintf "codes[%d] = window.atob(\"%s\");\n" i src);
    (* obtain smallstep diff *)
    if i > 0 then begin
      let diff = compute_diff() in
      out (sprintf "smallsteps[%d] = { exectime: %d, script: window.atob(\"%s\"), diff: window.atob(\"%s\") };\n" (i-1) stepdescr.exectime (Base64.encode_exn stepdescr.script) diff);
    end;
    (* obtain bigstep diff *)
    let isstartofbigstep =
      match stepdescr.isbigstep with
      | None -> false
      | Some _descr -> true
      in
    let isendofbigstep = (i > 0 && isstartofbigstep) || i = n-1 in
    if isendofbigstep then begin
      cmd "mv tmp_big.cpp tmp_before.cpp";
      let diff = compute_diff() in
      out (sprintf "bigsteps.push({ start: %d, stop: %d, descr: window.atob(\"%s\"), diff: window.atob(\"%s\") });\n"  !lastbigstepstart i (Base64.encode_exn !nextbigstep_descr) diff);
    end;
    (* shift files to prepare for next step *)
    cmd "mv tmp_after.cpp tmp_before.cpp";
    if isstartofbigstep then begin
      cmd "cp tmp_before.cpp tmp_big.cpp";
      lastbigstepstart := i;
      begin match stepdescr.isbigstep with
      | None -> assert false
      | Some descr -> nextbigstep_descr := descr
      end;
    end
  ) hist_and_descr;
  cmd "rm -f tmp.base64 tmp_after.cpp tmp_before.cpp tmp_big.cpp";
  close_out out_js
  *)


(* [dump_traces_to_js]: dump all traces to js.
    LATER: later generalize to multiple traces, currently it would
       probably overwrite the same file over and over again. *)
let dump_traces_to_js ?(prefix : string = "") () : unit =
  ()
  (* TODO FIX
  let history = get_decorated_history ~prefix () in
  dump_trace_to_js history
  *)

(*
  filename = prefix ^ "_trace.js"
  let f = open_out filename

  fprintf  f "var trace = {};";
  let print_step i ast =
     fprintf f "traces[%d] = {" i;
     code_to_js f i ast;
     fprintf f "};";
     in
  List.iteri print_step !traces



  var trace = {};
  trace[0] = {`
  trace[1] = {..};
  trace[2] = {..};
 *)


(* [step_tree_to_doc step_tree] takes a step tree and gives a string
   representation of it, using indentation to represent substeps *)
let step_tree_to_doc (step_tree:step_tree) : document =
  let ident_width = 3 in
  let rec aux (depth:int) (s:step_tree) : document =
    let i = s.step_infos in
    let space = blank 1 in
    let tab = blank (depth * ident_width) in
       tab
    ^^ string (step_kind_to_string s.step_kind)
    ^^ space
    ^^ string (sprintf "(%dms)" (int_of_float (i.step_duration *. 1000.)))
    ^^ space
    ^^ string i.step_name
    ^^ (separate space (List.map (fun (k,v) -> string (if k = "" then v else sprintf "~%s:%s" k v)) i.step_args))
    ^^ (if i.step_justif = [] then empty else separate empty (List.map (fun txt -> hardline ^^ tab ^^ string "==>" ^^ string txt) i.step_justif))
    ^^ (if i.step_script = "" then empty else (*hardline ^^ tab ^^*) string ">> " ^^ (string i.step_script))
    ^^ hardline
    ^^ separate empty (List.map (aux (depth+1)) s.step_sub)
    in
  aux 0 step_tree

(* [step_tree_to_file filename step_tree] takes a step tree and writes
   its string representation into a file. *)
let step_tree_to_file (filename:string) (step_tree:step_tree) =
  let line_width = 500 in
  let out = open_out filename in
  ToChannel.pretty 0.9 line_width out (step_tree_to_doc step_tree);
  close_out out

(* [dump_traces_to_textfile] dumps a trace into a text file *)
let dump_traces_to_textfile ?(prefix : string = "") () : unit =
  let prefix =
    if prefix = "" then the_trace.context.prefix else prefix in
  let filename = prefix ^ "_trace.txt" in
  printf "dumping trace to '%s'\n" filename;
  step_tree_to_file filename (get_root_step())


(* [output_prog_check_empty ~ast_and_enc ctx prefix ast_opt]: similar to [output_prog], but it
   generates an empty file in case the [ast] is an empty ast. *)
let output_prog_check_empty ?(ast_and_enc : bool = true) (ctx : context) (prefix : string) (ast_opt : trm) : unit =
  match ast_opt.desc with
  | Trm_seq tl when Mlist.length tl <> 0 -> output_prog ~ast_and_enc ctx prefix ast_opt
  | _ ->
      let file_prog = prefix ^ ctx.extension in
      let out_prog = open_out file_prog in
      close_out out_prog


(* [light_diff astBefore astAfter]: find all the functions that have not change after
    applying a transformation and hides their body for a more robust view diff. *)
let light_diff (astBefore : trm) (astAfter : trm) : trm * trm  =
    let topfun_before = top_level_fun_bindings astBefore in
    let topfun_after = top_level_fun_bindings astAfter in
    let topfun_common = get_common_top_fun topfun_before topfun_after in
    let filter_common ast = fst (hide_function_bodies (fun f -> List.mem f topfun_common) ast) in
    let new_astBefore = filter_common astBefore in
    let new_astAfter = filter_common astAfter in
    (new_astBefore, new_astAfter)

(* [dump_diff_and_exit()]: invokes [output_prog] on the current AST an also on the
   last item from the history, then it interrupts the execution of the script.
   This function is useful for interactively studying the effect of one particular
   small-step or big-step from the script. The diff is computed by taking the step
   at the top of the stack [the_trace], and considering the diff of the last sub-step
   the was performed. *)
(* LATER for mli: dump_diff_and_exit : unit -> unit *)
let dump_diff_and_exit () : unit =
  print_info None "Exiting script\n";
  let trace = the_trace in
  let ctx = trace.context in
  let prefix = (* ctx.directory ^ *) ctx.prefix in
  (* Common printinf function *)
  let output_ast ?(ast_and_enc:bool=true) filename_prefix ast =
    output_prog_check_empty ~ast_and_enc ctx filename_prefix ast;
    print_info None "Generated: %s%s\n" filename_prefix ctx.extension;
    in
  (* Extrat the two ASTs that should be used for the diff *)
  let step = get_cur_step() in
  let kind = step.step_kind in
  if kind <> Step_root && kind <> Step_big
    then failwith "dump_diff_and_exit: expects the current step to be a Root-step or a Big-step";
  let astBefore, astAfter =
    match step.step_sub with
    | [] -> failwith "dump_diff_and_exit: no sub-steps for which to display a diff";
    | last_step :: _ -> last_step.step_ast_before, last_step.step_ast_after
    in

  (* Option to compute light-diff:
      hide the bodies of functions that are identical in astBefore and astAfter. *)
  let astBefore, astAfter =
    if !Flags.use_light_diff
      then light_diff astBefore astAfter
      else astBefore, astAfter in

  (* Generate files *)
  output_ast (prefix ^ "_before") astBefore;
  output_ast (prefix ^ "_after") astAfter;
  print_info None "Writing ast and code into %s.js " prefix;
  (* Exit *)
  close_logs ();
  exit 0

(* [check_exit ~line] checks whether the program execution should be interrupted based
   on the command line argument [-exit-line]. If so, it exists the program after dumping
   the diff. *)
let check_exit ~(line:int) : unit (* may not return *) =
  let should_exit = match Flags.get_exit_line() with
    | Some li -> (line > li)
    | _ -> false
    in
  if should_exit
     then dump_diff_and_exit()

(* [check_exit_at_end] is called by [run.ml] when reaching the end of the script.
   This special case is needed to display a diff for the last transformation of
   the script. Indeed, this last transformation is not followed by a [!!] or a
   [bigstep] call.
   LATER: we may want to check that the targeted line is before the closing
   parenthesis of the call to script_cpp, but this requires instrumenting the
   call to script_cpp, and obtaining the end position of this construction. *)
let check_exit_at_end () : unit (* may not return *) =
  let should_exit = (Flags.get_exit_line () <> None) in
  if should_exit then begin
    close_smallstep_if_needed();
    if !Flags.only_big_steps
      then close_bigstep_if_needed();
    dump_diff_and_exit()
  end


(******************************************************************************)
(*                                   Steps                                     *)
(******************************************************************************)


(* [open_bigstep s]: announces that the next step is a bigstep, and registers
   a string description for that step. The [close_bigstep] is implicitly handled. *)
let open_bigstep ?(line : int = -1) (name:string) : unit =
  (* The [check_exit] is performed after closing the last small-step or big-step,
    depending on whether the user is interested in a diff over the last big-step. *)
  close_smallstep_if_needed();
  if not !Flags.only_big_steps
    then check_exit ~line;
  close_bigstep_if_needed();
  if !Flags.only_big_steps
    then check_exit ~line;
  (* Reparse if needed *)
  if !Flags.reparse_at_big_steps
    then reparse_alias ();
  open_step ~kind:Step_big ~name ~line ();
  (* Handle progress report *)
  if !Flags.report_big_steps then begin
    Printf.printf "Executing bigstep %s%s\n"
      (if line <> -1 then sprintf "at line %d" line else "")
      name
  end

(* [open_smallstep s]: announces that the next step is a smallstep,
   and registers a string description for that step, based on the excerpt
   frmo the file. The [close_smallstep] is implicitly handled. *)
(* LATER: add the line argument in the generation of the _with_lines file *)
let open_smallstep ?(line : int = -1) ?(reparse:bool=false) () : unit =
  close_smallstep_if_needed();
  if not !Flags.only_big_steps
    then check_exit ~line;
  if reparse
    then reparse_alias();
  let step_script =
    if !Flags.dump_trace
      then get_excerpt line
      else ""
    in
  open_step ~kind:Step_small ~name:"" ~line ~step_script ()


let transfo_step ~(name : string) ~(args : (string * string) list) (f : unit -> unit) : unit =
  step ~kind:Step_transfo ~name (fun () ->
    List.iter (fun (k, v) -> step_arg ~name:k ~value:v) args;
    f ()
  )

(* [check_recover_original()]: checks that the AST obtained so far
   is identical to the input AST, obtained from parsing. If not,
   it raises an error. *)
let check_recover_original () : unit =
  failwith "unimplemented"
  (*
  let check_same ast1 ast2 =
    if AstC_to_c.ast_to_string ast1 <> AstC_to_c.ast_to_string ast2
      then fail None "Trace.check_recover_original: the current AST is not identical to the original one."
      else () (* FOR DEBUG: Printf.printf "check_recover_original: successful" *)
    in
  let h = the_trace.history in
  match h with
  | [] -> failwith "check_recover_original: no history"
  | astLast :: [] -> () (* no operation performed, nothing to check *)
  | astLast :: astsBefore ->
      let _,astInit = Xlist.unlast astsBefore in
      check_same astLast astInit
*)



(******************************************************************************)
(*                                   User-level fucntions                     *)
(******************************************************************************)

  (* TODO: INTEGRATE Special hack for minimizing diff in documentation
  if !Flags.documentation_save_file_at_first_check <> "" then begin
    let trace = the_trace in
    let ctx = trace.context in
    output_prog ctx !Flags.documentation_save_file_at_first_check (trace.cur_ast)
 *)

(* [!!]: is a prefix notation for the operation [open_smallstep].
   By default, it performs only [step]. The preprocessor of the OCaml script file
   can add the [line] argument to the call to [open_smallstep], in order
   to allow for checking the exit line. Concretely, if the user has the cursor
   one line N when invoking the Optitrust "view_diff" command, then the tool
   will display the difference between the state of the AST at the first "!!"
   that occurs strictly after line N, and the state at the previous "!!",
   which could be on line N or before (or could correspond to the input AST
   loaded by [Trace.init] if there is no preceeding '!!'.).
   Use [!!();] for a step in front of another language construct, e.g., a let-binding. *)
let (!!) (x:'a) : 'a =
  open_smallstep ~reparse:false ();
  x

(* [!!!]: is similar to [!!] but forces a [reparse] prior to the [step] operation.
   ONLY FOR DEVELOPMENT PURPOSE. *)
let (!!!) (x : 'a) : 'a =
  open_smallstep ~reparse:true ();
  x


(* [dump ~prefix]: invokes [output_prog] to write the contents of the current AST.
   If there are several traces (e.g., due to a [switch]), it writes one file for each.
   If the prefix is not provided, the input file basename is used as prefix,
   and in any case "_out" is appended to the prefix.

   If you use [dump] in your script, make sure to call [!! Trace.dump] with the
   prefix [!!] in order for the diff visualization to work well for the last
   command before the call to dump.

   WILL BE DEPRECATED: If the command line argument [-dump-trace] was provided, then the
   function writes all the ASTs from the history into javascript files. *)
(* LATER for mli: val dump : ?prefix:string -> unit -> unit *)
let dump ?(prefix : string = "") () : unit =
  (* Dump final result, for every [switch] branch *)
  let ctx = the_trace.context in
  let prefix =
    if prefix = "" then (* ctx.directory ^ *) ctx.prefix else prefix
  in
  output_prog ctx (prefix ^ "_out") (the_trace.cur_ast)


(* DEPRECATED? [only_interactive_step line f]: invokes [f] only if the argument [line]
   matches the command line argument [-exit-line]. If so, it calls the
   [step] function to save the current AST, then calls [f] (for example
   to add decorators to the AST in the case of function [show]), then
   calls [dump_diff_and_exit] to visualize the effect of [f].
let only_interactive_step (line : int) ?(reparse : bool = false) (f : unit -> unit) : unit =
  let stepdescr_for_interactive_step =
    { isbigstep = None; script = ""; exectime = 0; } in
  if (Flags.get_exit_line() = Some line) then begin
    if reparse
      then
        reparse_alias ();
    step stepdescr_for_interactive_step;
    f();
    dump_diff_and_exit()
  end
  else
    begin
    check_exit_and_step();
    f()
    end
    *)

(* [ast()]: returns the current ast; this function should only be called within the
   scope of a call to [Trace.apply] or [Trace.call]. For example:
   [Trace.call (fun t -> ...  let t = ast() in ...) ].
   Note that in most cases, this function is not needed because the argument of
   the continuation already describes the current AST as the variable [t]. *)
let ast () : trm =
   the_trace.cur_ast

(* [set_ast]: is used for implementing [iteri_on_transformed_targets]. Don't use it elsewhere.
   NOTE: INTERNAL FUNCTION. *)
let set_ast (t:trm) : unit =
  the_trace.cur_ast <- t

(* [get_context ()]: returns the current context. Like [ast()], it should only be called
   within the scope of [Trace.apply] or [Trace.call]. *)
let get_context () : context =
  the_trace.context


(* LATER:  need to reparse to hide spurious parentheses *)
(* LATER: add a mechanism for automatic simplifications after every step *)
