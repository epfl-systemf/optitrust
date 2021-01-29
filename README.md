# Requirements

- OCaml (tested with version 4.08.0)
- dune (tested with version 2.5.1)
- [clangml](https://gitlab.inria.fr/tmartine/clangml) (tested with version
  4.2.0)
- clang-format (tested with version 6.0.0)
- [pprint](https://github.com/fpottier/pprint) (tested with version 20200410)
- Visual Studio Code (tested with version 1.44.1)
- Meld (tested with version 3.18.0)

# Installation details

```
   sudo apt-get install clang-format meld

   # Installation of opam: https://opam.ocaml.org/doc/Install.html
   sudo apt-get install opam
   opam switch create 4.11.0
   opam install dune clangml pprint

   # Installation of vscode: https://code.visualstudio.com/download
   # ...download the .deb package and install it

   # OCaml syntax highlighting for VSCode: https://marketplace.visualstudio.com/items?itemName=hackwaly.ocaml
   code &
   # Type CTRL+P, then paste and execute the commande:
   #    ext install hackwaly.ocaml
   # Follow tip1 from the link above to associate *.ml files with OCaml.
   # (optional) Disable minimap: menu "View" / uncheck "Show Minimap".

   # (optional) Install VSCode C++ extension
   #     ext install ms-vscode.cpptools
```

# Setup

The transformation script execution is based on a Visual Studio Code task. To
execute this task, one may use a keyboard shortcut.

Run `code` to open VSCode. To edit the `keybindings.json` file from Visual Studio
Code, type `Ctrl + Shift + P` to access the command panel and then choose
"Preferences: Open Keyboard Shortcuts (JSON)". There, replace the empty square 
braces with the following contents:

```
[
   {
       "key": "f6",
       "command": "workbench.action.tasks.runTask",
       "args": "Execute transformation script"
   },
   {
       "key": "ctrl+f6",
       "command": "workbench.action.tasks.runTask",
       "args": "Execute transformation script (dump trace)"
   }
]
```

This sets up `F6` to execute a transformation, and `CTRL+F6` to dump the full
trace of a transformation script (see `dump` in [SCRIPT.md](SCRIPT.md)). 


# Build and install

Execute `make && make install` at the root of the project.

# Example

```
  cp -R .vscode test_suite
  cd test_suite
  code . &
  # click on test_aosoa.ml in the list of files
  # press F6 to execute the comparison
```


# Usage

A transformation script is a `.ml` file. See [`SCRIPT.md`](SCRIPT.md) for
instructions on how to write scripts.

To edit a transformation script with the possibility of (partially) executing
it:
- Copy the `.vscode` directory into the directory that contains the script.
- In a terminal, go to the script directory and open Visual Studio Code by
  executing `code .`.
- Select the script file in the Visual Studio Code interface.
- To execute the script up to a given instruction, put the cursor on the line of
  this instruction and use the shortcut chosen during the setup phase.
  Assumption: the instruction holds on one line.
- To fully execute the script, just put the cursor on the line of the last
  instruction and use the shortcut.

When a script is (partially) executed, a diff for the source code before and
after the last transformation step is displayed with Meld.

# Source codes

Source codes in C or in C++ are allowed, but only a subset of these languages is
dealt with.

Constraints on switch statements:
- Cases must end with a break instruction.
- Nested cases must share their entire body. We call them case groups.

Constraint on `const` variables: it is forbidden to use the "address of"
operator on them.

Variables that are not `const` are heap allocated: the corresponding AST use
`const` pointers to such variables. This should be transparent for the user.

# Tests

Basic test scripts are provided in the following files:

- [`test_parser.ml`](test_suite/test_parser.ml).
- [`test_transformations.ml`](test_suite/test_transformations.ml).
- [`test_path.ml`](test_suite/test_path.ml).
- [`test_aosoa.ml`](test_suite/test_aosoa.ml).
- [`test_switch.ml`](test_suite/test_switch.ml).

To run tests, execute `make TESTS="test_name1 test_name2 ..." tests`, where
`test_name` is the name of a test script without its `.ml` extension. To run all
test scripts, simply execute `make tests`.
