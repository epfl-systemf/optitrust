open Optitrust
open Printf
module Terminal = Tools.Terminal

module StringSet = Set.Make(String)

(**

  This file is the source for the program `./tester.exe` used for executing unit tests.
  This program is meant to be called by the bash script `./tester`.
  It expects as arguments:

  1. the relative path from the optitrust root to the folder from which the user invoked the tester
  2. the name of an action to perform, e.g. 'run' (see further for other actions)
  3. other arguments, e.g. a list of test names, of test folder names.

  Usage: `./tester.exe [folder] [action] [options] [arg1] .. [argN]`.
  Usage:              `./tester [action] [options] [arg1] .. [argN]`.

  Options
  =======
  "-dry": only display list of tests to process
  "-hide-stdout": hide stdout contents produced by tests
  "-dump-trace": generate an html page describing every step performed for each test
  "-with-ignored": do not take into consideration the `ignored.tests`
  "-only-ignored": filter-out files that do not appear in an `ignored.tests` file
  "-yes": force the "yes" answer to actions requesting a confirmation
  "-no-doc": filter-out all files whose name contains the token "_doc."
  "-gitadd": automatically 'git add' the files generated by actions 'create' and 'addexp'
  "-no-lineshift": active this to help debug syntax errors in 'batch.ml'
  "-full-report": report a list telling for each test in order what was its status


  Action 'run'
  ==============
  Each [argi] must correspond to one of:
  - the full path to a `.ml` file, e.g., 'tests/loop/loop_fusion.ml', either relative to the root or relative to the caller's folder
  - a substring of a test name, e.g. 'fus', including:
    - the basename of a `.ml` file, e.g., 'loop_fusion.ml'
    - the basename without extension, e.g. 'loop_fusion'
  - the name of a folder appearing in 'tests/' or 'case_studies/', possibly in depth, e.g. 'loop'
   (LATER: should we restrict to subfolder of the caller's folder?)
  - the name of a '*.tests' file, e.g. 'fixme.tests'; particular case:
    if the name is 'ignore.tests', the option "-with-ignored" is automatically added

  If no [argi] is provided, then:
  - if [folder] is "." (the root of optitrust), then `all.tests` is used
    as the two [argi] arguments.
  - otherwise, [folder] is used as an [argi] argument.

  The tester program produces as output information on which tests fail.
  It generates a number of `*.tests` files:
  - `success.tests`: tests that have succeeded
  - `ignored.tests`: tests that have been ignored
  - `failed.tests`: tests whose execution raised an error
  - `missing_exp.tests`: tests not accompanied with a `_exp.cpp` file
  - `wrong.tests`: tests whose `_out.cpp` does not match the `_exp.cpp`.
  - `errors.tests`: tests that are not success (i.e. 'failed' or 'wrong' or 'missing_exp').
  These files are exploited by other actions meant for handling unsuccessful tests.

  In addition, the tester program updates the file `tofix.tests`:
  - if the file `tofix.tests` does not exist, it is initialized with the contents
    of `errors.tests`
  - if this file exists, it removes the lines that correspond tests that succeeded,
    and adding/merging the lines that correspond to tests that failed.

  The tester program ignores a `.ml` test file if there exists an `ignored.tests`
  file containing its name, in one of the folders containing the test file.
  The tests listed in the optional file 'wip_ignore.tests' (at the root)
  are also excluded by default.

  The optional file `wip_ignore.tests` is handled as followed:
  - tests in it are ignored like if they were in `ignore.tests`
    (unless the action is 'run+' instead of 'run')
  - for every test that is in `wip_ignore.tests` but not executed,
    a warning is printed at the end of the run, to remind the user
    that these tests still need to be fixed.

  In case a test is executed successfully, if it is mentioned in a 'ignore.tests'
  or 'wip_ignore.tests' file, then an info is printed to suggest the removal of
  the test from the ignore lists.


  Action 'create'
  ==============
  Each [argi] must be a path to a test to create, e.g. `tests/loop/loop_fusion.ml`.
  The `.ml` extension is optional.
  The tester program creates in that folder `loop_fusion.ml`, `loop_fusion.cpp`
  as well as `loop_fusion_doc.ml`, `loop_fusion_doc.cpp`.

  Action 'addexp'
  ==============
  If no [argi] is provided, the list of tests found in `missing_exp.tests` is used.
  For each test `foo` considered, the tester program checks that `foo_exp.cpp` does
  not yet exist, then copies the `foo_out.cpp` test to a `foo_exp.cpp`.

  Action 'fixexp'
  ==============
  If no [argi] is provided, the list of tests found in `wrong.tests` is used.
  For each test `foo` considered, the tester program copies checks that `foo_exp.cpp`
  exists, then copies the `foo_out.cpp` test to a `foo_exp.cpp`.

  Action 'ignore'
  ==============
  If no [argi] is provided, the list of tests found either in `failed.tests` or `wrong.tests` is used.
  For each test `foo` considered, the tester program adds the relative path to
  `foo.ml` as a line at the end of the `ignore.tests` located at the OptiTrust root.
  The purpose of this root `ignore.tests` is to act as a short-term buffer of broken
  tests, unlike the `ignore.tests` located in depth in the folders, which are meant
  list of tests that are meant to remain broken for a while.

  Action 'code'
  ==============
  If no [argi] is provided, the list of tests found in `failed.tests` or `wrong.tests` is used.
  For each test `foo` considered, the files `foo.ml` and `foo.cpp` are opened in VSCode.

  Action 'diff'
  ==============
  If no [argi] is provided, the list of tests found in `wrong.tests` is used.
  For each test `foo` considered, the command `code -d foo_out.cpp foo_exp.cpp`
  is executed.

  Action 'meld'
  ==============
  If no [argi] is provided, the list of tests found in `wrong.tests` is used.
  For each test `foo` considered, the command `meld foo_out.cpp foo_exp.cpp`
  is executed. Technically, a single call to `meld` is performed, opening
  all the pairs of files at once.

  ----
  LATER: an action 'sort':
  interactively navigate through the list of unsuccessful tests,
  (wrong.tests + failed.tests + missing_exp)
  and for each test ask an input key to decide wether to
  open in code, put in ignore, open diff, open meld, addexp, fixexp.

  LATER: an option to compare for each test in failed.tests the _out.cpp with the _exp.cpp,
  this is useful after modifying the _exp files by hand.

  TODO: run tofix.tests does not seem to exclude the ignored tests
*)

(*****************************************************************************)
(** Tools (LATER: move to tools/*.ml) *)

let debug = false

let debug_match_expected = ref false

exception TesterFailure of string

(* [fail msg] is for a fatal error produced by the tester *)
let fail (msg:string) : 'a =
  raise (TesterFailure msg)

let (~~) iter l f =
  iter f l


let do_is_ko (cmd : string) : bool =
  if debug then printf "%s\n" cmd;
  let exit_code = Sys.command cmd in
  exit_code != 0

let do_is_ok (cmd : string) : bool =
  not (do_is_ko cmd)

let do_or_die (cmd : string) : unit =
  if debug then printf "%s\n" cmd;
  let exit_code = Sys.command cmd in
  if exit_code != 0 then
    fail (sprintf "command '%s' failed with exit code '%i'" cmd exit_code)

  (* LATER: rediriiger des erreurs dans un fichier  2>&
    Sys.command en version booléenne
    ERRLINE = cat errorlog | head -n 1 | awk '{print $2}'
     ou grep ", line "
    head -n ${ERRLINE} batch.ml | grep "batching" | tail -1
     *)

(*****************************************************************************)
(** Options *)

(** Folders where tests might be located *)
(* LATER: remove and consider all mentioned directories instead *)
let tests_folders = ["tests"; "case_studies"]

(*** Folder from which the user invoked 'tester' *)
let caller_folder : string ref = ref ""

(* Action requested, e.g. 'run' *)
let action : string ref = ref ""

(** List of the [argi] arguments provided. *)
let args : string list ref = ref []

(** Flag to enable verbose mode *)
let verbose : bool ref = ref false

(** Flag to also generate traces *)
let dump_trace : bool ref = ref false

(** Flag to enable dry-run mode *)
let dry_run : bool ref = ref false

(** Flag to include tests appearing in ignored.tests files *)
let flags_with_ignored : bool ref = ref false

(** Flag to filter-out tests that do not appear in ignored.tests files *)
let flags_only_ignored : bool ref = ref false

(** Automatically git add new files *)
let git_add_new_files : bool ref = ref false

(** Flag to disable interactive confirmations *)
let skip_confirmation : bool ref = ref false

(** Flag to skip documentation-related tests *)
let skip_doc_tests : bool ref = ref false

(** Disable shifting of locations in batch.ml *)
let disable_lineshift : bool ref = ref false

(** Report the status of each test, in order *)
let full_report : bool ref = ref false


(*****************************************************************************)
(** FUTURE OPTIONS *) (* use of cache, and controlling output file generation *)

(* Flag for controlling whether or not to generate *_out.cpp files. *)
type outfile_gen =
  | Outfile_gen_always
  | Outfile_gen_only_on_failure (* failure or missing expected file *)
  | Outfile_gen_never

let outfile_gen : outfile_gen ref = ref Outfile_gen_only_on_failure

let string_to_outfile_gen = function
  | "always" -> Outfile_gen_always
  | "never" -> Outfile_gen_never
  | "onfailure" -> Outfile_gen_only_on_failure
  | _ -> fail "Invalid argument for -out"

let set_outfile_gen str =
  outfile_gen := string_to_outfile_gen str

(* Flag to control at which level the comparison is performed (AST or text).
   If Comparison_method_text, then implies Outfile_gen_always. *)
type comparison_method =
| Comparison_method_ast (* TODO: do we have this? *)
| Comparison_method_text

let comparison_method : comparison_method ref = ref Comparison_method_text

let _remove_later = comparison_method := Comparison_method_ast;
  comparison_method := Comparison_method_text

  (* Flag Comparison_method_text implies generation of output file
  if !comparison_method = Comparison_method_text
    then outfile_gen := Outfile_gen_always;
  *)

(*****************************************************************************)
(** Parsing of options *)

(* [cmdline_args]: a list of possible command line arguments. *)
type cmdline_args = (string * Arg.spec * string) list

(* [spec]: possible command line arguments. *)
let spec : cmdline_args =
   [ ("-dry", Arg.Set dry_run, " only display the list of tests to process");
     ("-hide-stdout", Arg.Set Flags.hide_stdout, " hide the contents that tests print on standard output ");
     ("-dump-trace", Arg.Set dump_trace, " generate html pages containing a trace of the steps performed by every test  ");
     ("-with-ignored", Arg.Set flags_with_ignored, " do not take into consideration the `ignore.tests` files  ");
     ("-only-ignored", Arg.Set flags_only_ignored, " filter-out files that do not appear in an `ignore.tests` file  ");
     ("-gitadd", Arg.Set git_add_new_files, " automatically call git add on files generated by `create` and `addexp`.");
     ("-v", Arg.Set verbose, " report details on the testing process.");
     ("-no-doc", Arg.Set skip_doc_tests, " skip tests whose names contains '_doc.' as a substring.");
     ("-yes", Arg.Set skip_confirmation, " answer 'yes' to confirmation prompts.");
     ("-no-lineshift", Arg.Set disable_lineshift, " disable shifting of lines in batch.ml.");
     ("-full-report", Arg.Set full_report, " report a list with the status of each test in order.");
     ("-all-warnings", Arg.Set Flags.report_all_warnings, " report all warnings.");
     ("-use-clang-format", Arg.Set Flags.use_clang_format, " forces the use of clang-format even in tester mode; automatically activated by -dump-trace.");

     (* NOT YET IMPLEMENTED *)
     ("-out", Arg.String set_outfile_gen, " generate output file: 'always', or 'never', or 'onfailure' (default)");

  ]


(*****************************************************************************)
(** Auxiliary function for handling dry-mode *)

(* [run_action] is like [do_or_die], with an option to print the command,
   and only prints the command if [-dry] flag has been set *)

let run_action ?(print = false) (cmd : string) : unit =
  (* FOR DEBUG:
   let _ = ignore print in
   printf "%s: %s\n" (if !dry_run then "[DRY]" else "[REAL]") cmd *)
  let pr () = printf "%s\n" cmd in
  if !dry_run then begin
    pr()
  end else begin
    if print then pr();
    do_or_die cmd
  end

(*****************************************************************************)
(** Processing of lists of tests provided as arguments *)

module File_set = Set.Make(String)

let file_set_ref_add r x =
  r := File_set.add x !r

let tmp_file = Filename.temp_file "command_output" ".txt"

let is_folder (file: string): bool =
  Sys.file_exists file && Sys.is_directory file

(** Return the list of tests mentioned in file [f], ignoring empty lines
    and lines commented out with a '#' as first character *)
let rec tests_from_file ?(relative = false) (file : string) : string list =
  if debug then printf "tests_from_file %s\n" file;
  let folder = Filename.dirname file in
  let lines = List.filter (fun s -> s <> "" && s.[0] <> '#') (Xfile.get_lines_or_empty file) in
  ~~ List.concat_map lines (fun test ->
    let test = if not relative || folder = "." then test else folder ^ "/" ^ test in
    (* TODO: this code should share common patterns with resolve_arg;
       there should almost no differences between an arg on the command line,
       and an arg in a .tests file *)
    if Filename.extension test = ".tests" then begin
      tests_from_file ~relative test
    end else if is_folder test then begin
      (* TODO: function for find all tests somewhere *)
      let sname = "-name '*.ml' -and -not -name '*_with_lines.ml'" in
      do_or_die (sprintf "find %s %s > %s" test sname tmp_file);
      Xfile.get_lines_or_empty tmp_file
    end else begin
      [test]
    end
  )

(** Call [f] on every test listed in file [f] *)
let foreach_test_from_file ?(relative = false) (file : string) (f : string -> unit) : unit =
  List.iter f (tests_from_file ~relative file)

(** Compute the list of all ignored tests listed in 'ignore.tests' files *)
let find_all_tests_to_ignore () : File_set.t =
  let target = String.concat " " tests_folders in
  (* LATER: use a version of Sys.command that captures the output *)
  do_or_die (sprintf "find %s -name 'ignore.tests' > %s" target tmp_file);
  let ignore_tests_files =
      (if Sys.file_exists "ignore.tests" then ["ignore.tests"] else [])
    @ (Xfile.get_lines_or_empty tmp_file) in
  if !verbose
    then printf "All ignore.tests files:\n  %s\n" (String.concat "\n  " ignore_tests_files);
  let result = ref File_set.empty in
  (* For each file named `ignore.tests` *)
  ~~ List.iter ignore_tests_files (fun ignore_tests_file ->
    (* For each test listed in that `ignore.tests`, ignoring lines starting with '#' *)
    foreach_test_from_file ~relative:true ignore_tests_file (fun ignore_test ->
        file_set_ref_add result ignore_test
    )
  );
  (* Return result *)
  if !verbose
    then printf "All tests ignored:\n  %s\n" (String.concat "\n  " (File_set.elements !result));
  !result

(** Resolve one [argi] to a list of tests (without filtering ignored tests),
    assuming its the name of a folder. *)
    (* LATER: do we want to deal with relative [folder]? *)
let resolve_arg_as_folder (arg : string) : string list =
  if is_folder arg then [arg] else begin
    let sfolders = String.concat " " tests_folders in
    do_or_die (sprintf "find %s -type d -name '%s' > %s" sfolders arg tmp_file);
    Xfile.get_lines_or_empty tmp_file
  end

(** [ensure_ml_extension path] applies to an absolute file path [path];
    if the filename suffix is [.cpp] or [_exp.cpp] or [_out.cpp], then it
    automatically changes the extension to the corresonding [.ml] file,
    and check that this files exist. *)
let ensure_ml_extension (path : string) : string =
  let error : 'a. unit -> 'a = fun () ->
    fail (sprintf "argument '%s' does not have the right extension." path) in
  let process (suffix : string) : string =
    let arg = (Tools.remove_suffix ~suffix path) ^ ".ml" in
    if not (Sys.file_exists arg)
      then error();
    arg in
  if String.ends_with ~suffix:".ml" path then
    path
  else if String.ends_with ~suffix:".cpp" path then
    process ".cpp"
  else if String.ends_with ~suffix:"_out.cpp" path then
    process "_out.cpp"
  else if String.ends_with ~suffix:"_exp.cpp" path then
    process "_exp.cpp"
  else
    error()

(** Resolve one [argi] to a list of tests (without filtering ignored tests) *)
let resolve_arg (arg : string) : string list =
  let is_file (fname : string) : bool =
    Sys.file_exists fname && not (Sys.is_directory fname) in
  let relarg = !caller_folder ^ "/" ^ arg in
  let extension = Filename.extension arg in
  if extension = ".tests" then begin
    if Filename.basename arg = "ignore.tests"
      then flags_with_ignored := true;
    if is_file relarg then begin
      if !verbose then printf "Resolved %s as relative .tests file: %s\n" arg relarg;
      tests_from_file ~relative:true relarg
    end else if is_file arg then begin
      if !verbose then printf "Resolved %s as absolute .tests file\n" arg;
      tests_from_file ~relative:true arg
    end else begin
      fail (sprintf "tester could not find file %s not %s" arg relarg)
    end
  end else if is_file relarg then begin
    (* Case is a relative file path *)
    if !verbose then printf "Resolved %s as relative filepath: %s\n" arg relarg;
    [ensure_ml_extension relarg]
  end else if is_file arg then begin
    (* Case an absolute file path *)
    if !verbose then printf "Resolved %s as absolute filepath\n" arg;
    [ensure_ml_extension arg]
  end else begin
    let (folders_to_search_from, pattern_on_the_name) : (string list * string option) =
      begin match resolve_arg_as_folder arg with
      | [] ->
        (* [arg] is not a folder, it is treated as a substring of a test name *)
        (* TODO: deal with relative [folder]. *)
        if !verbose then printf "Resolved %s as a substring\n" arg;
        tests_folders, Some arg
      | folders ->
        (* Case [arg] is exactly a directory name, possibly in depth *)
        if !verbose then printf "Resolved %s as folder name\n" arg;
        folders, None
      end in
    let sfolders = String.concat " " folders_to_search_from in
    let sname = "-name '*.ml' -and -not -name '*_with_lines.ml'" in
    let sname = match pattern_on_the_name with
      | None -> sname
      | Some pat -> sprintf "-name '*%s*' -and %s" pat sname
      in
    (* LATER: use a version of Sys.command that captures the output *)
    do_or_die (sprintf "find %s %s > %s" sfolders sname tmp_file);
    let tests = Xfile.get_lines_or_empty tmp_file in
    tests
  end

(** [check_test_extension test] checks that [test] has [.ml] has extension. *)
let check_test_extension (test : string) : unit =
  if Filename.extension test <> ".ml"
    then fail (sprintf "the test '%s' does not have .ml as extension" test)

(** [get_tests_and_ignored] Takes the [argi] arguments, and resolve them to tests,
   then filter out ignored tests based on the contents of all the `ignore.tests` files,
   each of which refers to tests with path relative to the location of the
   `ignore.tests` file that contains them.
   Influenced by the flags [-only-ignored] and [-with-ignored].
   Returns [tests_to_process, ignored_tests]. *)
let get_tests_and_ignored (args : string list) : (string list * string list) =
  let tests = File_set.of_list (List.concat_map resolve_arg args) in
  File_set.iter check_test_extension tests;
  let ignored = find_all_tests_to_ignore () in
  let tests_not_ignored = File_set.diff tests ignored in
  let tests_and_ignored = File_set.inter tests ignored in
  let set_to_process, set_to_ignore =
    if !flags_with_ignored then
      (tests, File_set.empty)
    else if !flags_only_ignored then
      (tests_and_ignored, tests_not_ignored)
    else
      (tests_not_ignored, tests_and_ignored) in
  let set_to_process, set_to_ignore =
    if not !skip_doc_tests then
      (set_to_process, set_to_ignore)
    else begin
      let is_doc (test : string) : bool =
        Tools.string_contains "_doc." (Filename.basename test) in
      let doc_tests, notdoc_tests = File_set.partition is_doc set_to_process in
      (notdoc_tests, File_set.union doc_tests set_to_ignore)
    end in
  File_set.elements set_to_process, File_set.elements set_to_ignore


(*****************************************************************************)
(** Action 'run' *)

(** Auxiliary function to compare the output with the expected output.
    Features an optimization to save the need for invoking clang-format,
    whenever possible. *)
let match_expected (filename_out:string) (filename_exp:string) : bool =
  if not (Sys.file_exists filename_out)
    then fail "match_expected expects filename_out to exit";
  if not (Sys.file_exists filename_exp)
    then fail "match_expected expects filename_exp to exit";
  let same_contents (f1 : string) (f2 : string) : bool =
    do_is_ok (sprintf "./tests/diff.sh %s %s > /dev/null" f1 f2) in
  if !Flags.use_clang_format then begin
    if !debug_match_expected then Tools.info "[use_clang_format] is on, cannot optimize [match_expected]";
    same_contents filename_out filename_exp
  end else begin
    (* At this point, [filename_out] is not formatted, whereas [filename_exp] is formatted;
       Both [orig_out] and [orig_exp], if they exist, are unformatted. *)
    let orig_out = Trace.filename_before_clang_format filename_out in
    let orig_exp = Trace.filename_before_clang_format filename_exp in
    let same = ref false in
    (* First, attempt to compare original files *)
    if Sys.file_exists orig_exp then begin
      if !debug_match_expected then Tools.info (sprintf "tested correctness without clang-format: %s" filename_out);
      let same_orig = same_contents filename_out orig_exp in
      if same_orig then same := true;
    end;
    (* If no original file, or if original files mismatch, compute clang-format and compare *)
    if not !same then begin
      (* format [filename_out], generating the file [orig_out] on the way *)
      Trace.cleanup_cpp_file_using_clang_format ~uncomment_pragma:true filename_out;
      if not (Sys.file_exists orig_out)
        then fail "match_expected assumes keep_file_before_clang_format to be true";
      (* now ready to compare formatted files *)
      let same_formatted = same_contents filename_out filename_exp in
      (* if debug_match_expected then Tools.info (sprintf "correctness is %s" (if same then "true" else "false"));*)
      if !debug_match_expected then Tools.info (sprintf "checked correctness using clang-format: %s" filename_out);
      if same_formatted then begin
        (* for next time, save the orig_out as orig_exp, because the result
          of formatting both files using clang-format is identical *)
        ignore (Sys.command (sprintf "cp %s %s" orig_out orig_exp));
        if !debug_match_expected then Tools.info (sprintf "saved unformatted output for future comparisons: %s." orig_exp);
      end;
      if same_formatted then same := true;
    end;
    !same
  end

(** Auxiliary function for updating the contents of 'tofix.tests',
    by removing contents from `success.tests`, or copying contents from
    'errors.tests' if 'tofix.tests' does not exist *)
let update_tofix_tests_list () : unit =
  let tofix_file = "tofix.tests" in
  let success_file = "success.tests" in
  let errors_file = "errors.tests" in
  if Sys.file_exists tofix_file then begin
    (* Remove successful *)
    do_or_die (sprintf "cat %s %s %s | sort | uniq --unique > %s"
      errors_file success_file success_file tofix_file);
    (* Add errors *)
    do_or_die (sprintf "cat %s >> %s; sort --unique %s -o %s"
      errors_file tofix_file tofix_file tofix_file)
  end else begin
    (* If nonexisting file, use contents of errors file *)
    do_or_die (sprintf "cp %s %s" errors_file tofix_file);
  end

let action_run (tests : string list) : unit =
  (* Compute tests to proces *)
  let (tests_to_process, tests_ignored) = get_tests_and_ignored tests in
  let nb_tests_to_process = List.length tests_to_process in
  let tests_to_process_string = String.concat " " tests_to_process in
  if !dry_run || !verbose
   then printf "Tester considering: \n  %s\nTester ignoring: \n  %s\nDry run completed.\n"
        (String.concat "\n  " tests_to_process)
        (String.concat "\n  " tests_ignored);
  if !dry_run then exit 0;
  if nb_tests_to_process = 0 then printf "Empty set of tests considered.\n";

  (* Always enable [-all-warnings] if there is only one test *)
  if nb_tests_to_process = 1
    then Flags.report_all_warnings := true;

  (* Enable backtrace display only when running an individual test *)
  if nb_tests_to_process > 1
    then Flags.print_backtrace_on_error := false;

  (* Generate a `batch.ml` program that contains the contatenation of the source
     code of every test considered *)
  (* LATER: could re-implement batch_tests.sh in OCaml *)
  let sdisable_lineshift = if !disable_lineshift then "export DISABLE_LINESHIFT=1; " else "" in
  do_or_die (sprintf "%stests/batch_tests.sh %s > tests/batch/batch.ml"
    sdisable_lineshift tests_to_process_string);

  (* Delete existing output files to avoid considering them in case an error occurs *)
  let delete_output test =
    let rm (f:string) : unit =
      (*if !debug_match_expected then Tools.info (sprintf "rm %s" f);*)
      ignore (do_is_ko (sprintf "rm -f %s > /dev/null" f)) in
    let prefix = Filename.remove_extension test in
    let filename_out = sprintf "%s_out.cpp" prefix in
    rm filename_out;
    let filename_orig_out = sprintf "%s_out_orig.cpp" prefix in
    rm filename_orig_out
    (* LATER: remove without displaying error messages or missing files *)
  in
  List.iter delete_output tests_to_process;

  (* Compile the `batch.ml` file, using dune hacks *)
  do_or_die "cp tests/batch/dune_disabled tests/batch/dune";
  let compile_success = do_is_ok "dune build tests/batch/batch.cmxs" in
  do_or_die "rm tests/batch/dune";
  if not compile_success then begin
    eprintf "failed to compile tests/batch/batch.ml";
    exit 2
  end;
  (* DEPRECATED printf "\n"; *)

  (* Start redirecting stdout into a temporary file during the execution
    of the unit tests. If -hide-stdout option is used, the contents of
    this file is ignored. Otherwise, the contents is printed at the end,
    or just before the error if an error is raised.  *)
  let oldstdout = Unix.dup Unix.stdout in
  let catpured_stdout_channel_ref = ref None in
  let catpured_stdout_file = Filename.temp_file "optitrust_batch_stdout" ".txt" in
  let capture = !Flags.hide_stdout in
  if capture then begin
    let catpured_stdout_channel = open_out catpured_stdout_file in
    Unix.dup2 (Unix.descr_of_out_channel catpured_stdout_channel) Unix.stdout;
    catpured_stdout_channel_ref := Some catpured_stdout_channel;
  end;
  let close_redirected_stdout () : unit =
    flush stdout;
    if capture then begin
      let catpured_stdout_channel = Option.get !catpured_stdout_channel_ref in
      close_out catpured_stdout_channel;
      Unix.dup2 oldstdout Unix.stdout;
      if not !Flags.hide_stdout then begin
        let catpured_stdout = Xfile.get_contents catpured_stdout_file in
        print_string catpured_stdout;
      end
    end in

  (* Execute the `batch.ml` program *)
  begin try
    Flags.program_name := "tester.ml";
    Dynlink.loadfile "tests/batch/batch.cmxs"
  with
  | Dynlink.Error err ->
      close_redirected_stdout();
      let sbt = Printexc.get_backtrace() in
      Printf.eprintf "%s\n%s" (Dynlink.error_message err) sbt;
      exit 3
  end;
  close_redirected_stdout();

  (* Analyse test results *)
  let tests_all = ref [] in
  let test_report (status : string) (acc : string list ref) (test : string) : unit =
    Tools.ref_list_add tests_all (sprintf "- %s\t%s" status test);
    Tools.ref_list_add acc test
    in
  let tests_failed = ref [] in
  let tests_noexp = ref [] in
  let tests_wrong = ref [] in
  let tests_success = ref [] in
  ~~ List.iter tests_to_process (fun test ->
    let prefix = Filename.chop_extension test in
    let filename_out = sprintf "%s_out.cpp" prefix in
    let filename_exp = sprintf "%s_exp.cpp" prefix in
    if Sys.file_exists filename_out then begin
      if Sys.file_exists filename_exp then begin
        (* TODO: use -q in the diff? *)
        if match_expected filename_out filename_exp
          then test_report "success" tests_success test
          else test_report "wrong" tests_wrong test
      end else begin
          test_report "noexp" tests_noexp test;
      end
    end else begin
      (* Missing _out.cpp file means test has failed *)
      test_report "failed" tests_failed test;
    end);

  (* Produce summary of errors, or full report of all tests *)
  let print_errors msg tests =
    if tests <> [] then
      printf "%s:\n%s\n" msg (Tools.list_to_string ~sep:"\n  " ~bounds:("  ", "\n") tests)
  in
  if !full_report then begin
    print_errors "All tests" !tests_all;
  end else begin
    print_errors "Failed tests" !tests_failed;
    print_errors "Missing expected" !tests_noexp;
    print_errors "Wrong tests" !tests_wrong;
  end;

  (* Produce .tests files *)
  Xfile.put_lines "success.tests" !tests_success;
  Xfile.put_lines "ignored.tests" tests_ignored;
  Xfile.put_lines "failed.tests" !tests_failed;
  Xfile.put_lines "wrong.tests" !tests_wrong;
  Xfile.put_lines "missing_exp.tests" !tests_noexp;
  Xfile.put_lines "errors.tests" (!tests_wrong @ !tests_failed @ !tests_noexp);
  update_tofix_tests_list();

  (* Produce general summary *)
  let isfirst = ref true in
  let print_count (color: string) (name: string) (tests: string list): unit =
    let len = List.length tests in
    if len > 0 then begin
      if not !isfirst then print_string ", ";
      print_string (Terminal.with_color color (sprintf "%i %s" len name));
      isfirst := false
    end
  in
  print_count Terminal.red "failed" !tests_failed;
  print_count Terminal.orange "missing exp" !tests_noexp;
  print_count Terminal.red "wrong" !tests_wrong;
  print_count Terminal.light_gray "ignored" tests_ignored;
  print_count Terminal.green "success" !tests_success;
  printf "\n";
  if !tests_failed <> [] || !tests_noexp <> [] || !tests_wrong <> [] then
    exit 1



(*****************************************************************************)
(** Action 'create' *)

let action_create (tests : string list) : unit =
  ~~ List.iter tests (fun test ->
    let prefix = Filename.remove_extension test in
    let extensions = [".ml"; ".cpp"; "_doc.ml"; "_doc.cpp"] in
    let files = List.map (fun s -> prefix ^ s) extensions in
    printf "Creating %s.{ml,cpp,_doc.ml,_doc.cpp}\n" prefix;
    ~~ List.iter files (fun file ->
      run_action (sprintf "touch %s" file);
      if !git_add_new_files then run_action (sprintf "git add %s" file);
      );
  )

(*****************************************************************************)
(** Action 'addexp' *)

let action_addexp (tests : string list) : unit =
  ~~ List.iter tests (fun test ->
    let prefix = Filename.remove_extension test in
    let outfile = prefix ^ "_out.cpp" in
    let expfile = prefix ^ "_exp.cpp" in
    (* TODO: the tester program checks that `foo_exp.cpp` does
  not yet exist *)
    run_action ~print:true (sprintf "cp %s %s" outfile expfile);
    if !git_add_new_files then run_action (sprintf "git add %s" expfile);
  )


(*****************************************************************************)
(** Action 'fixexp' *)

let action_fixexp (tests : string list) : unit =
  ~~ List.iter tests (fun test ->
    let prefix = Filename.remove_extension test in
    let outfile = prefix ^ "_out.cpp" in
    let expfile = prefix ^ "_exp.cpp" in
    if Sys.file_exists expfile
      then run_action ~print:true (sprintf "cp %s %s" outfile expfile)
      else eprintf "Warning: tester-fixexp does not find file %s\n" expfile;
  )


(*****************************************************************************)
(** Action 'ignore' *)

let action_ignore (tests : string list) : unit =
  ~~ List.iter tests (fun test ->
    let ignorefile = "./ignore.tests" in
    run_action (sprintf "echo '%s' >> %s" test ignorefile)
  )


(*****************************************************************************)
(** Action 'code' *)

let action_code (tests : string list) : unit =
  ~~ List.iter tests (fun test ->
    let prefix = Filename.remove_extension test in
    let mlfile = prefix ^ ".ml" in
    let cppfile = prefix ^ ".cpp" in
    run_action (sprintf "code %s %s" mlfile cppfile);
  )


(*****************************************************************************)
(** Action 'diff' *)

let action_diff (tests : string list) : unit =
  ~~ List.iter tests (fun test ->
    let prefix = Filename.remove_extension test in
    let outfile = prefix ^ "_out.cpp" in
    let expfile = prefix ^ "_exp.cpp" in
    run_action (sprintf "code -d %s %s" outfile expfile);
  )

(*****************************************************************************)
(** Action 'meld' *)

let action_meld (tests : string list) : unit =
  let meldargs = ref [] in
  ~~ List.iter tests (fun test ->
    let prefix = Filename.remove_extension test in
    let outfile = prefix ^ "_out.cpp" in
    let expfile = prefix ^ "_exp.cpp" in
    Tools.ref_list_add meldargs (sprintf "--diff %s %s " outfile expfile);
  );
  run_action (sprintf "meld %s" (String.concat "" !meldargs))


(*****************************************************************************)
(** Main *)

let _main : unit =
  (* Configure flags for tester *)
  (* [report_all_warnings] is false by default when using the tester on multiple tests;
     Can be changed with the [-all-warnings] option. *)
  Flags.report_all_warnings := false;
  (* [keep_file_before_clang_format] is set to true in order to allow for faster
     file comparisons. *)
  Flags.keep_file_before_clang_format := true;
  (* [use_clang_format] is false by default; unless producing traces *)
  let original_use_clang_format = !Flags.use_clang_format in
  Flags.use_clang_format := false;
  (* Trick to make sure that [ignore_serialized] is reset after each test to
     its original value. *)
  Flags.ignore_serialized_default := !Flags.ignore_serialized;

  (* Parsing of command line *)
  Arg.parse
    (Arg.align (spec @ Flags.spec))
    (fun arg ->
      if !caller_folder = ""
        then caller_folder := arg
      else if !action = ""
        then action := arg
      else
        Tools.ref_list_add args arg)
    "Usage: ./tester [action] [options] [arg1] .. [argN]\naction = run | create | addexp | fixexp | ignore | code | diff | meld\noptions:\n";

  (* Handle the dump_trace flag *)
  if !dump_trace then begin
    Flags.execution_mode := Execution_mode_full_trace;
    Flags.use_clang_format := original_use_clang_format;
  end;

  (* Check caller_folder has been provided *)
  if !caller_folder = ""
    then fail "a folder must be provided as first argument.";

  (* Check action has been provided *)
  if !action = ""
    then fail "an action must be provided as second argument.";

  (* Handle the case of an empty list of [argi], depending on [action] *)
  if !args = [] then begin
    if !action = "run" then begin
      if !caller_folder = "."
        then args := ["all.tests"]
        else args := [!caller_folder]
    end else if List.mem !action ["ignore"; "code"] then begin
      args :=   tests_from_file "failed.tests"
              @ tests_from_file "wrong.tests";
    end else if List.mem !action ["diff"; "meld"; "fixexp"] then begin
      args := tests_from_file "wrong.tests";
    end else if !action = "addexp" then begin
      args := tests_from_file "missing_exp.tests";
    end;
    if !verbose then printf "empty args resolved as: %s\n" (Tools.list_to_string !args)
  end;

  (* Switch according to action *)
  let execute () =
    let args = !args in
    match !action with
    | "run" -> action_run args
    | "create" -> action_create args
    | "addexp" -> action_addexp args
    | "fixexp" -> action_fixexp args
    | "ignore" -> action_ignore args
    | "code" -> action_code args
    | "diff" -> action_diff args
    | "meld" -> action_meld args
    | x -> fail (sprintf "unknown action %s" x)
    in

  (* Handle confirmation *)
  begin try
    let without_confirmation = ["run"; "meld"; "diff"] in
    let needs_confirmation = (not (List.mem !action without_confirmation)) && not !skip_confirmation in
    if not needs_confirmation then begin
      execute()
    end else begin
      let old_dry_run = !dry_run in
      dry_run := true;
      execute();
      dry_run := old_dry_run;
      if !dry_run then begin
        Tools.info "Dry run completed.";
        exit 0;
      end;
      flush stderr;
      printf "Confirm? Press ENTER to continue...\n";
      flush stdout;
      let answer = input_char stdin in
      if answer <> '\n' then begin
        Tools.info "Aborted.";
        exit 0;
      end;
      execute();
      Tools.info "Done."
    end
  with TesterFailure msg ->
    Tools.Terminal.(report red "TESTER ERROR" msg)
  end
