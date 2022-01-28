open Ast


(* [line_of_last_step] stores the line number from the source script at which a step
   ('!!' or '!^') was last processed. *)
let line_of_last_step = ref (-1)


(******************************************************************************)
(*                             Logging management                             *)
(******************************************************************************)

(* [timing_log] is a handle on the channel for writing timing reports. *)
let timing_log_handle = ref None

(* [logs] is a reference on the list of open log channels. *)
let logs : (out_channel list) ref = ref []

(* [close_logs] closes all open log channels. *)
let close_logs () : unit =
  List.iter (fun log -> close_out log) !logs;
  logs := []

(* [init_logs] initializes the log files. It closes any existing logs.
   Returns the one log created. *)
let init_logs directory prefix =
  close_logs();
  let clog = open_out (directory ^ prefix ^ ".log") in
  let timing_log = open_out ("timing.log") in
  timing_log_handle := Some timing_log;
  logs := timing_log :: clog :: [];
  clog

(* [write_log clog msg] writes the string [msg] to the channel [clog]. *)
let write_log (clog : out_channel) (msg : string) : unit =
  output_string clog msg;
  flush clog

(* [trm_to_log clog styp t] writes in the channel [clog] the term [t],
   and its typing information described by the string [styp]. *)
let trm_to_log (clog : out_channel) (exp_type : string) (t : trm) : unit =
  let sloc =
    match t.loc with
    | None -> ""
    | Some {loc_file = _; loc_start = {pos_line = start_row; pos_col = start_column}; loc_end = {pos_line = end_row; pos_col = end_column}} ->
       Printf.sprintf "at start_location %d  %d end location %d %d" start_row start_column end_row end_column
    in
  let msg = Printf.sprintf (" -expression\n%s\n" ^^ " %s is a %s\n") (Ast_to_c.ast_to_string t) sloc exp_type in
 write_log clog msg

(******************************************************************************)
(*                             Timing logs                                    *)
(******************************************************************************)

(* [write_timing_log msg] writes a message in the timing log file. *)
let write_timing_log (msg : string) : unit =
  let timing_log = match !timing_log_handle with
    | Some log -> log
    | None -> failwith "uninitialized timing log"
   in
   write_log timing_log msg

(* [timing ~name f] writes the execution time of [f] in the timing log file. *)
let timing ?(cond : bool = true) ?(name : string = "") (f : unit -> 'a) : 'a =
  if !Flags.analyse_time && cond then begin
    let t0 = Unix.gettimeofday () in
    let res = f() in
    let t1 = Unix.gettimeofday () in
    let msg = Printf.sprintf "%d\tms -- %s\n" (Tools.milliseconds_between t0 t1) name in
    write_timing_log msg;
    res
  end else begin
    f()
  end

(* [start_time] stores the date at which the script execution started (before parsing). *)
let start_time = ref (0.)

(* [last_time] stores the date at which the execution of the current step started. *)
let last_time = ref (0.)

(* [last_time_update()] returns the delay elapsed since the last call to this function;
   it updates the reference [last_time] with the current date. *)
let last_time_update () : int =
  let t0 = !last_time in
  let t = Unix.gettimeofday() in
  last_time := t;
  Tools.milliseconds_between t0 t

(* [report_time_of_last_step()] reports in the timing log the duration of the step
   that just completed, as measured by a call to [last_time_update]. *)
let report_time_of_last_step () : unit =
  if !Flags.analyse_time then begin
    let duration_of_previous_step = last_time_update () in
    write_timing_log (Printf.sprintf "===> TOTAL: %d\tms\n" duration_of_previous_step);
  end

(* [id_big_step] traces the number of big steps executed. This reference is used only
   when executing a script from the command line, because in this case the line numbers
   from the source script are not provided on calls to the [step] function. *)
let id_big_step = ref 0


(******************************************************************************)
(*                             File input                                     *)
(******************************************************************************)

(* [get_cpp_includes filename] get list of file includes syntactically visible
   on the first lines of a CPP file -- this implementation is quite restrictive. *)
let get_cpp_includes (filename : string) : string =
  (* make sure the include list is clean *)
  let includes = ref "" in
  let c_in = open_in filename in
  try
    while (true) do
      let s = input_line c_in in
      if Str.string_match (Str.regexp "^#include") s 0 then
        includes := !includes ^ s ^ "\n\n";
    done;
    !includes
  with
  | End_of_file -> close_in c_in; !includes

(* [parse filename] returns (1) a list of filenames corresponding to the '#include',
   and the OptiTrust AST. *)
let parse ?(parser = Parsers.Default) (filename : string) : string * trm =
  let use_new_encodings = !Flags.use_new_encodings in
  let parser = if parser = Parsers.Default then Flags.default_parser else parser in
  print_info None "Parsing %s...\n" filename;
  let includes = get_cpp_includes filename in
  let command_line_include =
    List.map Clang.Command_line.include_directory
      (Clang.default_include_directories ()) in
  let command_line_warnings = ["-Wno-parentheses-equality"; "-Wno-c++11-extensions"] in
  let command_line_args = command_line_warnings @ command_line_include in

  let t =
    timing ~name:"tr_ast" (fun () ->
      if use_new_encodings then begin
        let parse_clang () =
          Clang_to_astRawC.tr_ast (Clang.Ast.parse_file ~command_line_args filename) in
        let parse_menhir () =
          CMenhir_to_astRawC.tr_ast (MenhirC.parse_c_file_without_includes filename) in
        let rawAst = match parser with
          | Parsers.Default -> assert false (* see def of parser; Flags.default_parser should not be Default *)
          | Parsers.Clang -> parse_clang()
          | Parsers.Menhir -> parse_menhir()
          | Parsers.All ->
             let rawAstClang = parse_clang() in
             let rawAtMenhir = parse_menhir() in
             let strAstClang = Ast_to_rawC.ast_to_string rawAstClang in
             let strAstMenhir = Ast_to_rawC.ast_to_string rawAtMenhir in
             if strAstClang <> strAstMenhir then begin
               (* LATER: we could add a prefix based on the filename, but this is only for debug *)
               Xfile.put_contents "ast_clang.cpp" strAstClang;
               Xfile.put_contents "ast_menhir.cpp" strAstMenhir;
              fail None "parse: [-cparser all] option detected discrepencies; see ast_clang.cpp and ast_menhir.cpp";
             end else
             (* If the two ast match, we can use any one of them (only locations might differ); let's use the one from the default parser. *)
               if Flags.default_parser = Parsers.Clang then rawAstClang else rawAtMenhir
            in
          if !Flags.bypass_cfeatures
            then rawAst
            else CRawAst_to_ast.cfeatures_elim rawAst
      end else
        Clang_to_ast.translate_ast (Clang.Ast.parse_file ~command_line_args filename)
    )
  in
  (*  *)
  (* let ast =
    timing ~name:"parse_file" (fun () ->
      Clang.Ast.parse_file ~command_line_args filename
      ) in *)

  (* DEBUG: Format.eprintf "%a@."
       (Clang.Ast.format_diagnostics Clang.not_ignored_diagnostics) ast; *)
  print_info None "Parsing Done.\n";
  print_info None "Translating Done...\n";

  (* let t =
    timing ~name:"translate_ast" (fun () ->
      if use_new_encodings
        then Clang_to_astRawC.tr_ast ast
        else Clang_to_ast.translate_ast ast) in *)

  print_info None "Translation done.\n";
  (includes, t)

(******************************************************************************)
(*                             Trace management                               *)
(******************************************************************************)

(* A context contains general information about:
   - the source code that was loaded initially using [set_init_file],
   - the prefix of the filenames in which to output the final result using [dump]
   - the log file to report on the transformation performed. *)
type context =
  { extension : string;
    directory : string;
    prefix : string;
    includes : string;
    clog : out_channel; }

let context_dummy : context =
  { extension = ".cpp";
    directory = "";
    prefix = "";
    includes = "";
    clog = stdout; }

(* A trace is made of a context, a current AST, and a list of ASTs that were
   saved as "interesting intermediate steps", via the [Trace.save] function.
   Any call to the [step] function adds a copy of [cur_ast] into [history]. *)
type trace = {
  mutable context : context;
  mutable cur_ast : trm;
  mutable history : trms; }

let trm_dummy : trm =
  trm_val (Val_lit Lit_unit)

let trace_dummy : trace =
  { context = context_dummy;
    cur_ast = trm_dummy; (* dummy *)
    history = []; }

(* [traces] denotes the internal state of Optitrust. It consists of a list of traces,
   because Optitrust supports a [switch] command that allows branching in the
   transformation script, thus producing several possible traces. *)
type traces = trace list

let traces : traces ref =
  ref [trace_dummy]


(* [is_traces_dummy()] returns whether the trace was never initialized. *)
let is_traces_dummy () : bool =
  match !traces with
  | [tr] -> (tr == trace_dummy)
  | _ -> false

(* [reset()] restores the global state (object [traces]) in its uninitialized state,
   like at the start of the program. This operation is automatically called by [Trace.init]. *)
let reset () : unit =
  close_logs();
  traces := [trace_dummy]

(* [ml_file_excerpts] maps line numbers to the corresponding sections in-between [!!] marks in
   the source file. Line numbers are counted from 1 in that map. *)
module Int_map = Map.Make(Int)
let ml_file_excerpts = ref Int_map.empty

(* [compute_ml_file_excerpts lines] is a function for grouping lines according to the [!!] symbols. *)
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
  let regexp_double_quote = Str.regexp "^[ ]*!!" in
  let starts_with_double_quote (str : string) : bool =
    Str.string_match regexp_double_quote str 0 in
  let process_line (iline : int) (line : string) : unit =
    if starts_with_double_quote line then begin
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

let get_initial_ast (ser_mode : Flags.serialized_mode) (ser_file : string) (filename : string) : (string * trm) =
  (* LATER if ser_mode = Serialized_Make then let _ = Sys.command ("make " ^ ser_file) in (); *)
  let includes = get_cpp_includes filename in
  let ser_file_exists = Sys.file_exists ser_file in
  let ser_file_more_recent = if (not ser_file_exists) then false else Tools.is_file_newer_than ser_file filename in
  let auto_use_ser = (ser_mode = Serialized_Auto && ser_file_more_recent) in
  if (ser_mode = Serialized_Use
   || ser_mode = Serialized_Make
   || auto_use_ser) then begin
    if not ser_file_exists
      then fail None "get_initial_ast: please generate a serialized file first";
    if not ser_file_more_recent
      then fail None (Printf.sprintf "get_initial_ast: serialized file is out of date with respect to %s\n" filename);
    let ast = unserialize_from_file ser_file in
    if auto_use_ser
      then Printf.printf "Loaded ast from %s.\n" ser_file;
    (includes, ast)
    end
  else
    parse filename

(* [init f] initialize the trace with the contents of the file [f].
   This operation should be the first in a transformation script.
   The history is initialized with the initial AST.
   [~prefix:"foo"] allows to use a custom prefix for all output files,
   instead of the basename of [f]. *)
(* LATER for mli: val set_init_source : string -> unit *)
let init ?(prefix : string = "") (filename : string) : unit =
  reset ();
  let basename = Filename.basename filename in
  let extension = Filename.extension basename in
  let directory = (Filename.dirname filename) ^ "/" in
  let default_prefix = Filename.remove_extension basename in
  let ml_file_name =
    if Tools.pattern_matches "_inlined" default_prefix
      then List.nth (Str.split (Str.regexp "_inlined") default_prefix) 0
      else default_prefix in
  if !Flags.analyse_time then begin
    let src_file = (ml_file_name ^ ".ml") in
    if Sys.file_exists src_file then begin
      let lines = Xfile.get_lines src_file in
      ml_file_excerpts := compute_ml_file_excerpts lines;
    end;
  end;
  let mode = !Flags.serialized_mode in
  start_time := Unix.gettimeofday ();
  last_time := !start_time;
  let prefix = if prefix = "" then default_prefix else prefix in
  let clog = init_logs directory prefix in
  let ser_file = basename ^ ".ser" in
  let (includes, cur_ast) = get_initial_ast mode ser_file filename in
  let context = { extension; directory; prefix; includes; clog } in
  let trace = { context; cur_ast; history = [cur_ast] } in
  traces := [trace];
  if mode = Serialized_Build || mode = Serialized_Auto
    then serialize_to_file ser_file cur_ast;
  if mode = Serialized_Build
    then exit 0;
  print_info None "Starting script execution...\n"

(* [finalize()] should be called at the end of the script, to properly close the log files
    created by the call to [init]. *)
let finalize () : unit =
  close_logs()

(* [alternative f] executes the script [f] in the original state that
   was available just after the call to [init].
   After the call, all the actions performed are discarded.

  Current usage:
     !! Trace.alternative (fun () ->
        !! Loop.fusion_on_block [cLabel "tofusion"];
        !!());

   TODO: figure out if it is possible to avoid "!!" in front and tail of [Trace.restart].
   TODO: figure out if this implementation could be extended in the presence of [switch]. *)
let alternative f : unit =
  let saved_traces = !traces in
  let trace = match !traces with
    | [] -> fail None "alternative: the trace is empty"
    | [trace] -> trace
    | _ -> fail None "alternative: incompatible with the use of switch"
    in
  let init_ast =
    match List.rev trace.history with
    | [] -> fail None "alternative: the history is empty"
    | t::_ -> t
    in
  let init_trace = { trace with cur_ast = init_ast; history = [init_ast] } in
  traces := [init_trace];
  f();
  traces := saved_traces

(* [switch cases] allows to introduce a branching point in a script.
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
let switch ?(only_branch : int = 0) (cases : (unit -> unit) list) : unit =
  (* Close logs: new logs will be opened in every branch. *)
  close_logs ();
  let list_of_traces =
    Tools.fold_lefti
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

(* For dynamic checks: keep track of the number of nested calls to [Trace.call] *)
let call_depth = ref 0

(* [apply f] applies the transformation [f] to the current AST,
   and updates the current ast with the result of that transformation.
   If there are several active trace (e.g., after a [switch]),
   then [f] is applied to each of the traces. During the execution of [f]
   on a given trace, the set of traces is replaced with a singleton set
   made of only that trace; this allows for safe re-entrant calls
   (i.e., the function [f] itself may call [Trace.apply]. *)
let apply (f : trm -> trm) : unit =
  if is_traces_dummy()
    then fail None "Trace.init must be called prior to any transformation.";
  let cur_traces = !traces in
  incr call_depth;
  List.iter (fun trace ->
    traces := [trace]; (* temporary view on a single trace *)
    trace.cur_ast <- f trace.cur_ast)
    cur_traces;
  traces := cur_traces; (* restoring the original view on all traces *)
  decr call_depth

(* [call f] is similar to [apply] except that it applies to a function [f]
   with unit return type: [f] is meant to update the [cur_ast] by itself
   through calls to [apply].
   If there are several active trace (e.g., after a [switch]),
   then [f] is applied to each of the traces. During the execution of [f]
   on a given trace, the set of traces is replaced with a singleton set
   made of only that trace; this allows for safe re-entrant calls
   (i.e., the function [f] itself may call [Trace.apply]. *)
let call (f : trm -> unit) : unit =
  if is_traces_dummy()
    then fail None "Trace.init must be called prior to any transformation.";
  incr call_depth;
  let cur_traces = !traces in
  List.iter (fun trace ->
    traces := [trace]; (* temporary view on a single trace *)
    f trace.cur_ast)
    cur_traces;
  traces := cur_traces; (* restoring the original view on all traces *)
  decr call_depth

(* [step()] takes the current AST and adds it to the history.
   If there are several traces, it does so in every branch. *)
let step () : unit =
  List.iter (fun trace ->
    trace.history <- trace.cur_ast::trace.history)
    !traces

(* [check_recover_original()] checks that the AST obtained so far
   is identical to the input AST, obtained from parsing. If not,
   it raises an error. *)
let check_recover_original () : unit =
  let check_same ast1 ast2 =
    if Ast_to_rawC.ast_to_string ast1 <> Ast_to_rawC.ast_to_string ast2
      then fail None "check_recover_original: the current AST is not identical to the original one."
    in
  let check_trace trace =
    let h = trace.history in
    match h with
    | [] -> failwith "check_recover_original: no history"
    | astLast :: [] -> () (* no operation performed, nothing to check *)
    | astLast :: astsBefore ->
        let _,astInit = Tools.unlast astsBefore in
        check_same astLast astInit
    in
  List.iter check_trace !traces



(******************************************************************************)
(*                                   Output                                   *)
(******************************************************************************)

(* [cleanup_cpp_file_using_clang_format filename] makes a system call to
   reformat a CPP file using the clang format tool.
   LATER: find a way to remove extra parentheses in ast_to_doc, by using
   priorities to determine when parentheses are required. *)
let cleanup_cpp_file_using_clang_format (filename : string) : unit =
  timing ~name:(Printf.sprintf "cleanup_cpp_file_using_clang_format(%s)" filename) (fun () ->
    ignore (Sys.command ("clang-format -i " ^ filename)))

(* [output_prog ctx prefix ast] writes the program described by the term [ast]
   in several files:
   - one describing the raw AST ("prefix.ast")
   - one describing the internal AST ("prefix_enc.cpp")
   - one describing the CPP code ("prefix.cpp").
   The CPP code is automatically formatted using clang-format. *)


type language = | Language_cpp | Language_rust | Language_ocaml

let language_of_extension (extension:string) : language =
  match extension with
  | ".cpp" -> Language_cpp
  | ".rs" -> Language_rust
  | ".ml" -> Language_ocaml
  | _ -> fail None ("unknown extension " ^ extension)

let get_language () =
  match !traces with
  | [] -> fail None "cannot detect language -- trace should not be empty"
  | t::_ -> language_of_extension t.context.extension

let output_prog ?(beautify:bool=true) ?(ast_and_enc:bool=true) (ctx : context) (prefix : string) (ast : trm) : unit =
  let use_new_encodings = !Flags.use_new_encodings in
  let file_prog = prefix ^ ctx.extension in
  let out_prog = open_out file_prog in
  begin try
    (* print C++ code with decoding *)
    (*   DEPRECATED
    Printf.printf "===> %s \n" (ctx.includes); print_newline();*)
    (* LATER: try to find a way to put the includes in the AST so we can do simply ast_to_file *)
    output_string out_prog ctx.includes;
    if use_new_encodings then begin
      if !Flags.bypass_cfeatures
        then Ast_to_rawC.ast_to_outchannel ~optitrust_syntax:true out_prog ast
        else Ast_to_rawC.ast_to_outchannel out_prog (CRawAst_to_ast.cfeatures_intro ast)
    end else
      Ast_to_c.ast_to_outchannel out_prog ast;
    output_string out_prog "\n";
    close_out out_prog;
  with | Failure s ->
    close_out out_prog;
    failwith s
  end;
  (* beautify the C++ code --comment out for debug *)
  if beautify
    then cleanup_cpp_file_using_clang_format file_prog;
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
      output_string out_enc ctx.includes;
      if use_new_encodings then ()
        (* then Ast_to_rawC.ast_to_outchannel out_prog ast *)
        else Ast_to_c.ast_to_undecoded_doc out_enc ast;
      output_string out_enc "\n";
      close_out out_enc;
      cleanup_cpp_file_using_clang_format file_enc;
    with | Failure s ->
      close_out out_ast;
      close_out out_enc;
      failwith s
    end
  end

(* [output_prog_opt ?ast_and_enc ctx prefix ast_opt] is similar to [output_prog], but it
   generates an empty file in case the [ast_opt] is [None]. *)
let output_prog_opt ?(ast_and_enc : bool = true) (ctx : context) (prefix : string) (ast_opt : trm option) : unit =
  match ast_opt with
  | Some ast -> output_prog ~ast_and_enc ctx prefix ast
  | None ->
      let file_prog = prefix ^ ctx.extension in
      let out_prog = open_out file_prog in
      close_out out_prog

(* [output_js index cpp_filename prefix _ast] Produces a javascript file used for viewing the ast in interactive mode
   This javascript file contains an array of source codes and an array of ast's. Where the entry at index i contains the state
   of the source and ast after applying transformaion i.
*)
let output_js  ?(vars_declared : bool = false)(index : int) (prefix : string) (ast : trm) : unit =
  let file_js = prefix ^ ".js" in
  let out_js = open_out file_js in
  try
    (* Dump the description of the AST nodes *)
    let lang, extension = match get_language () with
      | Language_cpp -> "\'" ^ "text/x-c++src" ^ "\'", ".cpp"
      | Language_rust -> "\'" ^ "text/x-rustsrc" ^ "\'", ".rs"
      | Language_ocaml -> "\'" ^ "text/x-Ocaml" ^ "\'", ".ml" in
    if not vars_declared
      then begin
      output_string out_js (Tools.sprintf "var source = %s\n" "new Array();");
      output_string out_js (Tools.sprintf "var contents = %s\n" "new Array();");
      output_string out_js (Tools.sprintf "var language = %s\n" lang);
      end else ();
    let src = Xfile.get_contents (prefix ^ "_before" ^ extension) in
    Ast_to_js.Json.code_to_js out_js index src;
    output_string out_js "\n";
    Ast_to_js.ast_to_js out_js index ast;
    output_string out_js "\n";
    close_out out_js;
  with | Failure s ->
    close_out out_js;
    failwith s

(* [dump_trace_to_js] writes into one/several (?) files
   the contents of the current AST and of all the history,
   that is, of all the ASTs for which the [step] method was called.
   DEPRECATED, will be fixed soon. *)
let dump_trace_to_js ?(prefix : string = "") () : unit =
  assert (prefix = prefix && false);
  let dump_history (prefix : string) (asts : trms) : unit =
    let nbAst = List.length asts in
    let i = ref (nbAst - 2) in
    List.iter
      (fun ast ->
        if !i = 0 then
          output_js 0 prefix ast
        else if !i = (nbAst - 2) then
          output_js ~vars_declared:true (nbAst - 2) prefix ast
        else if !i = (-1) then ()
        else
          begin
          output_js ~vars_declared:true !i prefix ast;
          i := !i - 1
          end
      )
      asts in
  List.iter
    (fun trace ->
      let ctx = trace.context in
      let prefix =
        if prefix = "" then ctx.directory ^ ctx.prefix else prefix
      in
      dump_history prefix (trace.cur_ast :: trace.history)
    )
    (!traces)

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

(******************************************************************************)
(*                                   Reparse                                  *)
(******************************************************************************)

(* [reparse_trm ctx ast] print [ast] in a temporary file and reparses it using Clang. *)
let reparse_trm ?(info : string = "") ?(parser = Parsers.Default) (ctx : context) (ast : trm) : trm =
  if !Flags.debug_reparse then begin
    let info = if info <> "" then info else "of a term during the step starting at" in
    Printf.printf "Reparse: %s line %d.\n" info !line_of_last_step;
    flush stdout
  end;
  let in_prefix = ctx.directory ^ "tmp_" ^ ctx.prefix in
  output_prog ~beautify:false ctx in_prefix ast;

  let (_, t) = parse ~parser (in_prefix ^ ctx.extension) in
  (*let _ = Sys.command ("rm " ^ in_prefix ^ "*") in*)
  t

(* [reparse()] function takes the current AST, prints it to a file, and parses it
   as if it was a fresh input. Doing so ensures in particular that all the type
   information is properly set up. WARNING: reparsing discards all the marks in the AST. *)
let reparse ?(info : string = "") ?(parser = Parsers.Default) () : unit =
 List.iter (fun trace ->
    let info = if info <> "" then info else "the code during the step starting at" in
    trace.cur_ast <- reparse_trm ~info ~parser trace.context trace.cur_ast)
    !traces

(* Work-around for a name clash *)
let reparse_alias = reparse


(******************************************************************************)
(*                                   Dump                                     *)
(******************************************************************************)

(* [light_diff astBefore astAfter] finds all the functions that have not change after
    applying a transformation and hides their body for a more robust view diff. *)
 (* TODO: see the comment in dump_diff_and_exit on how to remove the option to simplify the code *)
let light_diff (astBefore : trm option) (astAfter : trm) : trm option * trm  =
  match astBefore with
  | Some astBefore ->
    let topfun_before = top_level_fun_bindings astBefore in
    let topfun_after = top_level_fun_bindings astAfter in
    let topfun_common = get_common_top_fun topfun_before topfun_after in
    let filter_common ast = fst (hide_function_bodies (fun f -> List.mem f topfun_common) ast) in
    let new_astBefore = filter_common astBefore in
    let new_astAfter = filter_common astAfter in
    (Some new_astBefore, new_astAfter)
  | _ -> astBefore, astAfter

(* [dump_diff_and_exit()] invokes [output_prog] on the current AST an also on the
   last item from the history, then it interrupts the execution of the script.
   This function is useful for interactively studying the effect of one particular
   transformation from the script.
   If option [-dump-last nb] was provided, output files are produced for the last [nb] step. *)
(* LATER for mli: dump_diff_and_exit : unit -> unit *)
let dump_diff_and_exit () : unit =
  if !Flags.analyse_time then begin
     report_time_of_last_step();
     write_timing_log (Printf.sprintf "------------------------TOTAL TRANSFO TIME: %.3f s\n" (!last_time -. !start_time));
     write_timing_log (Printf.sprintf "------------START DUMP------------\n");
  end;
  timing ~name:"TOTAL for dump_diff_and_exit" (fun () ->
    print_info None "Exiting script\n";
    let trace =
      match !traces with
      | [] -> fail None "No trace"
      | [tr] -> tr
      | trs -> Printf.eprintf "Warning: considering the last branch of all switches.\n";
              List.hd (List.rev trs)
      in
    let ctx = trace.context in
    let prefix = ctx.directory ^ ctx.prefix in
    (* Common printinf function *)
    let output_ast ?(ast_and_enc:bool=true) filename_prefix ast_opt =
      output_prog_opt ~ast_and_enc ctx filename_prefix ast_opt;
      print_info None "Generated: %s%s\n" filename_prefix ctx.extension;
      in
    (* CPP and AST output for BEFORE *)
    (* TODO: we could simplify quite a bit all this code by creating an empty AST
       when trace.history is empty. In ast.ml, you can define [empty_ast] to be
       a trm_seq with no items, and with the annotation Main_file. *)
    let astBefore =
      match trace.history with
      | t::_ -> Some t (* the most recently saved AST *)
      | [] -> Printf.eprintf "Warning: only one step in the history; consider previous step blank.\n"; None
      in
    let astAfter = trace.cur_ast in

    (* Compute light-diff: hide bodies of functions that are identical in astBefore and astAfter. *)
    let astBefore, astAfter =
      if !Flags.use_light_diff then light_diff astBefore astAfter else astBefore, astAfter in

    output_ast (prefix ^ "_before") astBefore;
    (* CPP and AST for BEFORE_N *)
    if !Flags.dump_last <> Flags.dump_last_default then begin
      let nb_requested = !Flags.dump_last in
      let nb_available = List.length trace.history in
      (* if nb_requested < nb_available
        then Printf.eprintf "Warning: not enought many steps for [dump_last]; completing with blank files.\n"; *)
      for i = 0 to nb_requested-1 do
        let astBeforeI = if i < nb_available then Some (List.nth trace.history i) else None in
        output_ast ~ast_and_enc:false (prefix ^ "_before_" ^ string_of_int i) astBeforeI
      done;
    end;
    (* CPP and AST and Javscript for AFTER *)
    output_ast (prefix ^ "_after") (Some astAfter);
    print_info None "Writing ast and code into %s.js " prefix;
    output_js 0 prefix astAfter;
    (* Printf.printf "EXIT   %s\n" prefix; *)
  );
  (* Exit *)
  close_logs ();
  exit 0

(* [check_exit_and_step()] performs a call to [check_exit], to check whether
   the program execution should be interrupted based on the command line argument
   [-exit-line], then it performas a call to [step], to save the current AST
   in the history, allowing for a visualizing the diff if the next call to
   [check_exit_and_step] triggers a call to [dump_diff_and_exit].
   If the optional argument [~reparse:true] is passed to the function,
   then the [reparse] function is called, replacing the current AST with
   a freshly parsed and typechecked version of it.
   The [~is_small_step] flag indicates whether the current step is small
   and should be ignored when visualizing big steps only. *)
let check_exit_and_step ?(line : int = -1) ?(is_small_step : bool = true)  ?(reparse : bool = false) () : unit =
  (* Update the line of the last step entered *)
  line_of_last_step := line;
  (* Special hack for minimizing diff in documentation *)
  if !Flags.documentation_save_file_at_first_check <> "" then begin
    let trace =
      match !traces with
      | [trace] -> trace
      | _ -> fail None "doc_script_cpp: does not support the use of [switch]"
      in
    let ctx = trace.context in
    output_prog ctx !Flags.documentation_save_file_at_first_check (trace.cur_ast)
  end else begin
    let ignore_step = is_small_step && !Flags.only_big_steps in
    if not ignore_step then begin
      report_time_of_last_step();
      (* Handle exit of script *)
      let should_exit =
        match Flags.get_exit_line() with
        | Some li -> (line > li)
        | _ -> false
        in
      if should_exit then begin
        if !Flags.analyse_time then begin
          write_timing_log (Printf.sprintf "------------------------\n");
        end;
        dump_diff_and_exit();
      end else begin
        (* Handle reparse of code *)
        if reparse || (!Flags.reparse_at_big_steps && not is_small_step) then begin
          let info = if reparse then "the code on demand at" else "the code just before the big step at" in
          let parser = !Flags.use_parser in
          reparse_alias ~info ~parser ();
          if !Flags.analyse_time then
            let duration_of_reparse = last_time_update () in
            write_timing_log (Printf.sprintf "------------------------\nREPARSE: %d\tms\n" duration_of_reparse);
        end;
        (* Handle the reporting of the execution time *)
        if !Flags.analyse_time then begin
          let txt =
            if !ml_file_excerpts = Int_map.empty then "" else begin
              match Int_map.find_opt line !ml_file_excerpts with
              | Some txt -> txt
              | None -> (*failwith*) Printf.sprintf "<unable to retrieve line %d from script>" line
            end in
          write_timing_log (Printf.sprintf "------------------------\n[line %d]\n%s\n" line txt);
        end;
        (* Handle progress report *)
        if not is_small_step && !Flags.report_big_steps then begin
          if line = -1 then begin
            incr id_big_step;
            Printf.printf "Executing big-step #%d\n" !id_big_step
          end else begin
            Printf.printf "Executing big-step line %d\n" line
          end;
          flush stdout
        end;
      end;
      (* Save the current code in the trace *)
      step();
  end
end

(* [!!] is a prefix notation for the operation [check_exit_and_step].
   By default, it performs only [step]. The preprocessor of the OCaml script file
   can add the [line] argument to the call to [check_exit_and_step], in order
   to allow for checking the exit line. Concretely, if the user has the cursor
   one line N when invoking the Optitrust "view_diff" command, then the tool
   will display the difference between the state of the AST at the first "!!"
   that occurs strictly after line N, and the state at the previous "!!",
   which could be on line N or before (or could correspond to the input AST
   loaded by [Trace.init] if there is no preceeding '!!'.).
   Use [!!();] for a step in front of another language construct, e.g., a let-binding. *)
let (!!) (x:'a) : 'a =
  check_exit_and_step ~is_small_step:true ~reparse:false ();
  x

(* [!^] is similar to [!!] but indicates the start of a big step in the transformation script. *)
let (!^) (x:'a) : 'a =
  check_exit_and_step ~is_small_step:false ~reparse:false ();
  x

(* [!!!] is similar to [!!] but forces a [reparse] prior to the [step] operation.
   ONLY FOR DEVELOPMENT PURPOSE. *)
let (!!!) (x : 'a) : 'a =
  check_exit_and_step ~is_small_step:true ~reparse:true ();
  x

(* [!!^] is forces reparse before a big step.
   ONLY FOR DEVELOPMENT PURPOSE. *)
let (!!^) (x : 'a) : 'a =
  check_exit_and_step ~is_small_step:false ~reparse:true ();
  x

(* [dump ~prefix] invokes [output_prog] to write the contents of the current AST.
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
  if !Flags.analyse_time then begin
      write_timing_log (Printf.sprintf "------------START DUMP------------\n");
  end;
  (* Dump full trace if requested *)
  if !Flags.dump_all then
     dump_trace_to_js ~prefix (); (* dump_trace ~prefix () *)
  (* Dump final result, for every [switch] branch *)
  List.iter
    (fun trace ->
      let ctx = trace.context in
      let prefix =
        if prefix = "" then ctx.directory ^ ctx.prefix else prefix
      in
      output_prog ctx (prefix ^ "_out") (trace.cur_ast)
    )
    (!traces)

(* [only_interactive_step line f] invokes [f] only if the argument [line]
   matches the command line argument [-exit-line]. If so, it calls the
   [step] function to save the current AST, then calls [f] (for example
   to add decorators to the AST in the case of function [show]), then
   calls [dump_diff_and_exit] to visualize the effect of [f]. *)
let only_interactive_step (line : int) ?(reparse : bool = false) (f : unit -> unit) : unit =
  if (Flags.get_exit_line() = Some line) then begin
    if reparse
      then
        let parser = !Flags.use_parser in
        reparse_alias ~parser ();
    step();
    f();
    dump_diff_and_exit()
  end
  else
    begin
    check_exit_and_step();
    f()
    end

(* [ast()] returns the current ast; this function should only be called within the
   scope of a call to [Trace.apply] or [Trace.call]. For example:
   [Trace.call (fun t -> ...  let t = ast() in ...) ].
   Note that in most cases, this function is not needed because the argument of
   the continuation already describes the current AST as the variable [t]. *)
let ast () : trm =
  if !call_depth = 0
    then failwith "[get_the_ast] can only be invoked inside a call to [Trace.call].";
   match !traces with
   | [tr] -> tr.cur_ast
   | [] -> assert false (* [!traces] can never be empty *)
   | _ -> failwith "[get_the_ast] can only be invoked inside a call to [Trace.call] and not after a switch."

(* INTERNAL FUNCTION.
   [set_ast] is used for implementing [iteri_on_transformed_targets]. Don't use it elsewhere. *)
let set_ast (t:trm) : unit =
  assert (!call_depth > 0);
  match !traces with
  | [tr] -> tr.cur_ast <- t
  | _ -> assert false

(* [get_context ()] returns the current context. Like [ast()], it should only be called
   within the scope of [Trace.apply] or [Trace.call]. *)
let get_context () : context =
  match !traces with
  | [tr] -> tr.context
  | _ -> fail None "get_context: couldn't get the current context"

(* [parse_cstring] is a function that can be used to acquire the AST of a statement
   provided as a string by the user, e.g, [int x; x = 1; x = 5;]. The arguments are:
   - [ctx]: an optional context, from which to obtain the list of '#include' to use.
   - [context]: describes a context in which the statement can be parsed
     (e.g., to parse [int x = y], in a context where [int y] is defined.
   - [is_expression]: a flag to indicate if we are parsing an expression or a statement
     (for expressions, we add a semicolon at the end).
   - [s]: the string that describes the code of which we want the AST. *)
(* TODO: change [ctx : context]   to [ctx option]. If it is Some, includes [ctx.includes]
   in addition  to the string [context]. *)
let parse_cstring (context : string) (is_expression : bool) (s : string) (ctx : context) : trms =
 let context = if context = "" then ctx.includes else context in
 let command_line_args =
  List.map Clang.Command_line.include_directory
    (ctx.directory :: Clang.default_include_directories())
  in
 let ast =
    Clang.Ast.parse_string ~command_line_args
      (Printf.sprintf
         {|
          %s
          void f(void){
            #pragma clang diagnostic ignored "-Wunused-value"
            %s
          }
          |}
         context
         (if is_expression then s ^ ";" else s)
      )
  in
  let t = Clang_to_ast.translate_ast ast in
  match t.desc with
  | Trm_seq tl1 when Mlist.length tl1 = 1 ->
    let t = Mlist.nth tl1 0 in
     begin match t.desc with
     | Trm_seq tl  ->
        let fun_def = List.nth (List.rev (Mlist.to_list tl)) 0 in
        begin match fun_def.desc with
        | Trm_let_fun (_, _, _, fun_body) ->
          begin match fun_body.desc with
          | Trm_seq tl -> Mlist.to_list tl
          | _ -> fail fun_body.loc "parse_cstring: expcted a sequence of terms"
          end
        | _ -> fail fun_def.loc "parse_cstring: expected a function definition"
        end
     | _ -> fail t.loc "parse_cstring: expected another sequence"
     end
  | _-> fail t.loc (Printf.sprintf "parse_cstring: exptected a sequence with only one trm, got %s\n" (Ast_to_c.ast_to_string t))

(* [Trace.term ctx s] returns the AST that corresponds to a statement
  described by the string [s]. The context in which the statement is
  parsed can be provided as optional argument. By default, we use only
  as context the '#include' obtained from the [ctx], which may be obtained
  using [get_context()]. *)
  (* TODO: should make ctx on optional argument ?(ctx:context), which could be None
    or Some ctx, in which case we would use the include from the context. *)
let term ?(context : string = "") (ctx : context) (s : string) : trm =
  let tl = parse_cstring context true s ctx  in
  match tl with
  | [expr] -> expr
  | _ -> fail None "term: expcted a list with only one element"

(* LATER:  need to reparse to hide spurious parentheses *)
(* LATER: add a mechanism for automatic simplifications after every step *)


