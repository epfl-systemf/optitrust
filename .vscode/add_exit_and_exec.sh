#!/bin/bash
# Arguments:
#   0. file dirname
#   1. filename of the transformation script (without extension)
#   2. active line number (to add exit_script instruction)
#   3. option(s) for execution
#      currently: only -dump-trace

# first step: add exit_script instruction at the end of the active line

# NOTE: if vscode is invoked from the test_suite folder, then no need for test_suite

DIRNAME=$1
FILEBASE=$2
LINE=$3
OPTIONS=$4

cd ${DIRNAME}
ocaml .vscode/add_exit.ml -file "${FILEBASE}.ml" -line ${LINE}

# second step: build and execute the script
ocamlbuild -pkgs clangml,refl,pprint,str,optiTrust.scriptTools "${FILEBASE}_with_exit.byte" || (echo "Cannot compile $1_with_exit.ml"; exit 1)
./${FILEBASE}_with_exit.byte ${OPTIONS}

# third step: clean up and show the diff of the two last states of the program
ocamlbuild -clean
rm "${FILEBASE}_with_exit.ml"
