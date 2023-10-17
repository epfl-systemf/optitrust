open Optitrust
open Printf

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
  "-no-doc": filter-out all files whose name contains the token "_doc_"


  Action 'run'
  ==============
  Each [argi] must correspond to one of:
  - the full path to a `.ml` file, e.g., 'tests/loop/loop_fusion.ml', either relative to the root or relative to the caller's folder
  - a substring of a test name, e.g. 'fus', including:
    - the basename of a `.ml` file, e.g., 'loop_fusion.ml'
    - the basename without extension, e.g. 'loop_fusion'
  - the name of a folder appearing in 'tests/' or 'case_studies/', possibly in depth, e.g. 'loop'
   (LATER: should we restrict to subfolder of the caller's folder?)
  - the name 'ignore.tests', in this case the option "-with-ignored" is automatically added
  - the name of a '*.tests' file, e.g. 'fixme.tests';



  If no [argi] is provided, then:
  - if [folder] is "." (the root of optitrust), then `tests case_studies` are used
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

  In addition, the tester program updates the file `tofix.tests`, if this file exists, by removing from it the lines that correspond to a test listed in `success.tests`. If the file `tofix.tests` does not exist, it is initialized with the contents of `errors.tests`.

  The tester program ignores a `.ml` test file if there exists an `ignored.tests` file containing its name, in one of the folders containing the test file.

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
*)

(*****************************************************************************)
(** Tools (LATER: move to tools/*.ml) *)

let ref_list_add (r: 'a list ref) (x: 'a) : unit =
  r := x :: !r

let (~~) iter l f =
  iter f l

let do_is_ko (cmd : string) : bool =
  let exit_code = Sys.command cmd in
  exit_code != 0

let _do_is_ok (cmd : string) : bool =
  not (do_is_ko cmd)

let do_or_die (cmd : string) : unit =
  let exit_code = Sys.command cmd in
  if exit_code != 0 then
    failwith (sprintf "command '%s' failed with exit code '%i'" cmd exit_code)

  (* LATER: rediriiger des erreurs dans un fichier  2>&
    Sys.command en version booléenne
    ERRLINE = cat errorlog | head -n 1 | awk '{print $2}'
     ou grep ", line "
    head -n ${ERRLINE} batch.ml | grep "batching" | tail -1
     *)

(*****************************************************************************)
(** Options *)

(** Folders where tests might be located *)
let tests_folders = ["tests"; "case_studies"]

(*** Folder from which the user invoked 'tester' *)
let caller_folder : string ref = ref ""

(* Action requested, e.g. 'run' *)
let action : string ref = ref ""

(** List of the [argi] arguments provided. *)
let args : string list ref = ref []

(** Flag to enable verbose mode *)
let verbose : bool ref = ref false

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
  | _ -> failwith "Invalid argument for -out"

let set_outfile_gen str =
  outfile_gen := string_to_outfile_gen str

(* Flag to ignore all cached data *)
let ignore_cache : bool ref = ref false

(* Flag to discard all cached data *)
let discard_cache : bool ref = ref false

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
     ("-dump-trace", Arg.Set Flags.dump_trace, " generate html pages containing a full trace of the steps performed by every test  ");
     ("-with-ignored", Arg.Set flags_with_ignored, " do not take into consideration the `ignore.tests` files  ");
     ("-only-ignored", Arg.Set flags_only_ignored, " filter-out files that do not appear in an `ignore.tests` file  ");
     ("-gitadd", Arg.Set git_add_new_files, " automatically call git add on files generated by `create` and `addexp`.");
     ("-v", Arg.Set verbose, " report details on the testing process.");
     (* NOT YET IMPLEMENTED *)
     ("-out", Arg.String set_outfile_gen, " generate output file: 'always', or 'never', or 'onfailure' (default)");
     ("-ignore-cache", Arg.Set ignore_cache, " ignore the serialized AST, force reparse of source files; does not modify the existing serialized data");
     ("-discard-cache", Arg.Set discard_cache, " clear all serialized AST; save serizalize data for tests that are executed.");
     ("-no-doc", Arg.Set skip_doc_tests, " skip tests whose names contains '_doc_' as a substring.");
     ("-yes", Arg.Set skip_confirmation, " answer 'yes' to confirmation prompts.");
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

let tmpfile = Filename.temp_file "command_output" ".txt"

(** Return the list of tests mentioned in file [f], ignoring empty lines *)
let tests_from_file (file : string) : string list =
  List.filter (fun s -> s <> "") (Xfile.get_lines_or_empty file)

(** Call [f] on every test listed in file [f] *)
let foreach_test_from_file (file : string) (f : string -> unit) : unit =
  List.iter f (tests_from_file file)

(** Compute the list of all ignored tests listed in 'ignore.tests' files *)
let find_all_tests_to_ignore () : File_set.t =
  let target = String.concat " " tests_folders in
  (* LATER: use a version of Sys.command that captures the output *)
  do_or_die (sprintf "find %s -name 'ignore.tests' > %s" target tmpfile);
  let ignore_tests_files =
      (if Sys.file_exists "ignore.tests" then ["ignore.tests"] else [])
    @ (Xfile.get_lines_or_empty tmpfile) in
  if !verbose
    then eprintf "All ignore.tests files:\n  %s\n" (String.concat "\n  " ignore_tests_files);
  let result = ref File_set.empty in
  (* For each file named `ignore.tests` *)
  ~~ List.iter ignore_tests_files (fun ignore_tests_file ->
    let folder = Filename.dirname ignore_tests_file in
    (* For each test listed in that `ignore.tests`, ignoring lines starting with '#' *)
    foreach_test_from_file ignore_tests_file (fun ignore_test ->
      if ignore_test <> "" && ignore_test.[0] <> '#' then begin
        let test = if folder = "." then ignore_test else folder ^ "/" ^ ignore_test in
        file_set_ref_add result test
      end
    )
  );
  (* Return result *)
  if !verbose
    then eprintf "All tests ignored:\n  %s\n" (String.concat "\n  " (File_set.elements !result));
  !result

(** Resolve one [argi] to a list of tests (without filtering ignored tests),
    assuming its the name of a folder. *)
    (* LATER: do we want to deal with relative [folder]? *)
let resolve_arg_as_folder (arg : string) : string list =
  if Sys.file_exists arg && Sys.is_directory arg then [arg] else begin
    let sfolders = String.concat " " tests_folders in
    do_or_die (sprintf "find %s -type d -name '%s' > %s" sfolders arg tmpfile);
    Xfile.get_lines_or_empty tmpfile
  end

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
      if !verbose then eprintf "Resolved %s as relative .tests file: %s\n" arg relarg;
      tests_from_file relarg
    end else if is_file arg then begin
      if !verbose then eprintf "Resolved %s as absolute .tests file\n" arg;
      tests_from_file arg
    end else begin
      failwith (sprintf "tester could not find file %s not %s" arg relarg)
    end
  end else if is_file relarg then begin
    (* Case is a relative file path *)
    if !verbose then eprintf "Resolved %s as relative filepath: %s\n" arg relarg;
    [relarg]
  end else if is_file arg then begin
    (* Case an absolute file path *)
    if !verbose then eprintf "Resolved %s as absolute filepath\n" arg;
    [arg]
  end else begin
    let (folders_to_search_from, pattern_on_the_name) : (string list * string option) =
      begin match resolve_arg_as_folder arg with
      | [] ->
        (* [arg] is not a folder, it is treated as a substring of a test name *)
        (* TODO: deal with relative [folder]. *)
        if !verbose then eprintf "Resolved %s as a substring\n" arg;
        tests_folders, Some arg
      | folders ->
        (* Case [arg] is exactly a directory name, possibly in depth *)
        if !verbose then eprintf "Resolved %s as folder name\n" arg;
        folders, None
      end in
    let sfolders = String.concat " " folders_to_search_from in
    let sname = "-name '*.ml' -and -not -name '*_with_lines.ml'" in
    let sname = match pattern_on_the_name with
      | None -> sname
      | Some pat -> sprintf "-name '*%s*' -and %s" pat sname
      in
    (* LATER: use a version of Sys.command that captures the output *)
    do_or_die (sprintf "find %s %s > %s" sfolders sname tmpfile);
    let tests = Xfile.get_lines_or_empty tmpfile in
    tests
  end

(** Takes the [argi] arguments, and resolve them to tests, then filter out
   ignored tests based on the contents of all the `ignore.tests` files,
   each of which refers to tests with path relative to the location of the
   `ignore.tests` file that contains them.
   Influenced by the flags [-only-ignored] and [-with-ignored].
   Returns [tests_to_process, ignored_tests]. *)
let get_tests_and_ignored (args : string list) : (string list * string list) =
  if !skip_doc_tests then failwith "option -no-doc is not yet supported"; (* TODO *)
  let tests = File_set.of_list (List.concat_map resolve_arg args) in
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
  File_set.elements set_to_process, File_set.elements set_to_ignore


(*****************************************************************************)
(** Action 'run' *)

(** Auxiliary function for updating the contents of 'tofix.tests',
    by removing contents from `success.tests`, or copying contents from
    'errors.tests' if 'tofix.tests' does not exist *)
let update_tofix_tests_list () : unit =
  let tofix_file = "tofix.tests" in
  let success_file = "success.tests" in
  let errors_file = "errors.tests" in
  if Sys.file_exists tofix_file then begin
    do_or_die (sprintf "cat %s %s %s | sort | uniq --unique > %s"
      errors_file success_file success_file tofix_file)
  end else begin
    do_or_die (sprintf "cp %s %s" errors_file tofix_file);
  end

let action_run (tests : string list) : unit =
  (* Compute tests to proces *)
  let (tests_to_process, tests_ignored) = get_tests_and_ignored tests in
  let nb_tests_to_process = List.length tests_to_process in
  let tests_to_process_string = String.concat " " tests_to_process in
  if !dry_run || !verbose
   then eprintf "Tester considering: \n  %s\nTester ignoring: \n  %s\nDry run completed.\n"
        (String.concat "\n  " tests_to_process)
        (String.concat "\n  " tests_ignored);
  if !dry_run then exit 0;
  if nb_tests_to_process = 0 then eprintf "Empty set of tests considered.";

  (* Enable backtrace display only when running an individual test *)
  if nb_tests_to_process > 1
    then Flags.print_backtrace_on_error := false;

  (* Generate a `batch.ml` program that contains the contatenation of the source
     code of every test considered *)
  (* LATER: could re-implement batch_tests.sh in OCaml *)
  do_or_die (sprintf "tests/batch_tests.sh %s > tests/batch/batch.ml" tests_to_process_string);

  (* Delete existing output files to avoid considering them in case an error occurs *)
  let delete_output test =
    let prefix = Filename.remove_extension test in
    let filename_out = sprintf "%s_out.cpp" prefix in
    ignore (do_is_ko (sprintf "rm -f %s > /dev/null" filename_out))
    (* LATER: remove without displaying error messages or missing files *)
  in
  List.iter delete_output tests_to_process;

  (* Compile the `batch.ml` file, using dune hacks *)
  do_or_die "cp tests/batch/dune_disabled tests/batch/dune";
  do_or_die "dune build tests/batch/batch.cmxs; rm tests/batch/dune";
  (* DEPRECATED printf "\n"; *)

  (* If -hide-stdout option is used, start redirecting stdout into
     "_tests_stdout.txt" during the execution of the unit tests *)
  let oldstdout = Unix.dup Unix.stdout in
  let newstdout = ref None in
  if !Flags.hide_stdout then begin
    let c = open_out "_tests_stdout.txt" in
    Unix.dup2 (Unix.descr_of_out_channel c) Unix.stdout;
    newstdout := Some c;
  end;
  let close_redirected_stdout () : unit =
    if !Flags.hide_stdout then begin
      flush stdout;
      let c = Option.get !newstdout in
      close_out c;
      Unix.dup2 oldstdout Unix.stdout;
    end in

  (* Execute the `batch.ml` program *)
  begin try
    Flags.program_name := "tester.ml";
    Dynlink.loadfile "tests/batch/batch.cmxs"
  with
    Dynlink.Error err -> begin
      close_redirected_stdout();
      let sbt = Printexc.get_backtrace() in
      Printf.eprintf "%s\n%s" (Dynlink.error_message err) sbt;
      exit 1
    end
  end;
  close_redirected_stdout();

  (* Analyse test results *)
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
        let mismatches_expected = do_is_ko (sprintf "./tests/diff.sh %s %s > /dev/null" filename_out filename_exp) in
        if mismatches_expected
          then ref_list_add tests_wrong test
          else ref_list_add tests_success test
      end else begin
          ref_list_add tests_noexp test;
      end
    end else begin
      (* Missing _out.cpp file means test has failed *)
      ref_list_add tests_failed test;
    end);

  (* Produce summary of errors *)
  let print_errors msg tests =
    if tests <> [] then
      printf "%s:\n%s\n" msg (Tools.list_to_string ~sep:"\n  " ~bounds:[""; "\n"] tests)
  in
  print_errors "Failed tests" !tests_failed;
  print_errors "Missing expected" !tests_noexp;
  print_errors "Wrong tests" !tests_wrong;

  (* Produce .tests files *)
  Xfile.put_lines "success.tests" !tests_success;
  Xfile.put_lines "ignored.tests" tests_ignored;
  Xfile.put_lines "failed.tests" !tests_failed;
  Xfile.put_lines "wrong.tests" !tests_wrong;
  Xfile.put_lines "missing_exp.tests" !tests_noexp;
  Xfile.put_lines "errors.tests" (!tests_wrong @ !tests_failed @ !tests_noexp);
  update_tofix_tests_list();

  (* Produce general summary *)
  let print_count (name: string) (tests: string list): unit =
    let len = List.length tests in
    if len > 0 then printf "%i %s, " len name
  in
  print_count "failed" !tests_failed;
  print_count "missing exp" !tests_noexp;
  print_count "wrong" !tests_wrong;
  print_count "ignored" tests_ignored;
  print_count "success" !tests_success;
  printf "\n"


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
    ref_list_add meldargs (sprintf "--diff %s %s " outfile expfile);
  );
  run_action (sprintf "meld %s" (String.concat "" !meldargs))


(*****************************************************************************)
(** Main *)

let _main : unit =
  (* Parsing of command line *)
  Arg.parse
    (Arg.align (spec @ Flags.spec))
    (fun arg ->
      if !caller_folder = ""
        then caller_folder := arg
      else if !action = ""
        then action := arg
      else
        ref_list_add args arg)
    "Usage: TODO";

  (* Check caller_folder has been provided *)
  if !caller_folder = ""
    then failwith "Invalid usage: a folder must be provided as first argument.";

  (* Check action has been provided *)
  if !action = ""
    then failwith "Invalid usage: an action must be provided as second argument.";

  (* Handle the case of an empty list of [argi], depending on [action] *)
  if !args = [] then begin
    if !action = "run" then begin
      if !caller_folder = "."
        then args := tests_folders
        else args := [!caller_folder]
    end else if List.mem !action ["ignore"; "code"] then begin
      args :=   tests_from_file "failed.tests"
              @ tests_from_file "wrong.tests";
    end else if List.mem !action ["diff"; "meld"; "fixexp"] then begin
      args := tests_from_file "wrong.tests";
    end else if !action = "addexp" then begin
      args := tests_from_file "missing_exp.tests";
    end
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
    | x -> failwith (sprintf "Invalid usage: unknown action %s" x)
    in

  (* Handle confirmation *)
  let needs_confirmation = (!action <> "run") && not !skip_confirmation in
  if not needs_confirmation then begin
    execute()
  end else begin
    let old_dry_run = !dry_run in
    dry_run := true;
    execute();
    dry_run := old_dry_run;
    if !dry_run then begin
      printf "Dry run completed.\n";
      exit 0;
    end;
    flush stderr;
    printf "Confirm? Press ENTER to continue...\n";
    flush stdout;
    let answer = input_char stdin in
    if answer <> '\n' then begin
      printf "Aborted.\n";
      exit 0;
    end;
    execute();
    printf "Done.\n"
  end


(*****************************************************************************)
(* DEPREACTED/FUTURE *)
  (*
     Produire une liste de (testname, result)

     type result =
       | Result_failure of string * string    --> short and full descr
       | Result_match
       | Result_mismatch

     if only one failure, print full descr for this one;
     in more, generate a file with all full_descr

     print all succeeded first
       print all mismatch (short failure descr => one per line)
       print all failed execs
     at last, print a summary: ALL SUCCEED, NB OF EXEC/DIFF FAILED, NB OF SKIPPED

     generate a file "failed.txt" with the list of tests that have Result_mismatch,
     for use with bash script for accepting changes

     LATER: generate a file "report.txt" or both with the list of tests that have succeeded

     *)

     (* NOTE:

Xfile.
let serialize_to_file (filename : string) (obj : 'a) : unit =
let unserialize_from_file (filename : string) : 'a =
let is_newer_than (filename1 : string) (filename2 : string) : bool =

need to compare dependency on
- optitrust/ast.cmxa installed needs to be more recent than serialized files
*)

(* Level 2 :
   if optitrust lib
     AND test.cpp
     AND test_out.cpp
     AND test.ml
     have not changed since production of report.ser
   then skip this test.
*)

(* --- accept.sh

  for each test in failed.txt
    meld test_out.cpp test_exp.cpp
    if diff test_out.cpp test_exp.cpp remains nonempty
    echo "accept change ? Y / N / A"
    x = input
    si Y, faire un cp
    si N, continue
    si A, exit

*)

(* -
NOT YET IMPLEMENTED
- caching of AST representation for input and expected output files
- management of dependencies on the source code of optitrust
- selection of a subset of tests to process (${args} is by default "all")
  (either via filenames, or via a 'key' name refering to a groupe of files)
- comparison either at the AST level or at the cpp source level via diff
----------
*)


  (* TODO: We cache the "raw ast".
  from a trm t, we need to serialize Ast_fromto_AstC.cfeatures_intro t;
  this is what is provided in trace.ml to  the function
   AstC_to_c.ast_to_outchannel ~beautify_mindex ~comment_pragma:use_clang_format out_prog t_after_cfeatures_intro;

    LATER: deal with script_cpp ~filename by searching 'batch.ml'

  for each test:
    test.cpp must exit
    if test.ser exist and .ser up to date: OK
    else update .ser by parsing it from .cpp and serialize

  for each test:
    if test_exp.cpp exist:
      if test_exp.ser exist and  .ser up to date: OK
      else update .ser by parsing it from _exp.cpp and serialize
    else:
      create test_exp.ser with empty ast
  *)
  (* Need to save the cached_inputs and cached_expoutputs :
     -> we want save the ones that were serialized before,
        and update the ones that have just been reparsed *)

(* c'est le code de batch.ml
     qui fait la gestion des cached_inputs/cached_outputs

     et qui pourrait faire la comparaison au niveau AST

     ./batch_controller.ml prendrait en argument presque tous les arguments de tester.ml

     tester.ml would do only the computation of the list of tests
     and the generation of batch.ml and the launching of batcher.exe

  *)


(*****************************************************************************)
(** Processing of tests list DEPRECATED *)
(*

let filename_concat folder filename =
  if folder = "." then filename else Filename.concat folder filename

let get_alias_targets (alias_filename: string) : string list =
  let folder = Filename.dirname alias_filename in
  let lines = Xfile.get_lines_or_empty alias_filename in
  List.filter_map (fun l ->
      if String.starts_with ~prefix:"#" l || String.length l = 0
        then None
        else Some (filename_concat folder l)) lines

let get_tests_in_dir (folder: string) : string list * File_set.t =
  let ignored_files = get_alias_targets (filename_concat folder "ignored.tests") in
  let ignored_files = List.fold_left (fun acc f -> File_set.add f acc) File_set.empty ignored_files in
  let folder_files = Sys.readdir folder in
  let test_files = List.filter_map (fun f ->
      if String.ends_with ~suffix:".ml" f && not (String.ends_with ~suffix:"_with_lines.ml" f)
      then
        let filename = filename_concat folder f in
        if File_set.mem filename ignored_files
          then None
          else Some filename
      else None) (Array.to_list folder_files) in
  List.sort String.compare test_files, ignored_files

let rec resolve_test_targets (target_list: string list) : string list * File_set.t =
  let test_files, ignored_file_set =
    List.fold_right (fun target (test_files, ignored_files) ->
    if String.ends_with ~suffix:".ml" target then
      (target :: test_files, ignored_files)
    else
      let alias_filename = target ^ ".tests" in
      let (new_test_files, new_ignored_files) =
        if Sys.file_exists alias_filename then
          let sub_targets = get_alias_targets alias_filename in
          resolve_test_targets sub_targets
        else
          get_tests_in_dir target
      in
      (new_test_files @ test_files, File_set.union ignored_files new_ignored_files)
  ) target_list ([], File_set.empty)
  in
  (test_files, ignored_file_set)


(* Takes the list of target arguments on the command line;
   and expand the 'targets', remove duplicates and ignored tests.
   Returns (tests_to_process, ignored_tests)
*)
let get_tests_and_ignored (targets: string list): (string list * string list) =
  let test_files, ignored_file_set = resolve_test_targets targets in
  let test_files = Xlist.remove_duplicates test_files in
  (* Tests that were otherwise selected are not ignored *)
  let ignored_file_set =
    List.fold_left (fun acc t -> File_set.remove t acc) ignored_file_set test_files
  in
  let _check_all_files_exist =
    List.iter (fun test_file ->
      if not (Sys.file_exists test_file)
        then failwith (sprintf "File not found:"))
      test_files;
    in
  (test_files, File_set.elements ignored_file_set)
*)


(*

      RES=$(find tests -name "*${arg}*.ml" -and -not -name "*_with_lines.ml")
*)
