#!/bin/bash


# TODO FOR ARTHUR: take as argument FILENAME, then deduce FILEBASE AND FILEEXT
# AND if we're calling a cpp file, then don't add_exit

DIRNAME=$1
FILEBASE=$2
LINE=$3
VIEW=$4 # should be view_diff or view_result
RECOMPILE_OPTITRUST=$5 # should be recompile_optitrust_yes or recompile_optitrust_no
OPTIONS=$6
OPTIONS2=$7

# LATER: if a ${FILEBASE}_exp.cpp file is present, export it into the JS file,
# so that the browser can report on the differences between _out.cpp and _exp.cpp.

#UPDATE=noupdate

# Path to .vscode folder and src folder and src/src folder
VSCODE=`pwd`
SRCFOLDER=`cd .. && pwd`
SRCSRCFOLDER=`cd ../src && pwd`

# This can help with opam switches
eval $(opam env)

# Make sure we work in the directory that contains the file
cd ${DIRNAME}


# Read options from optitrust_flags.sh
# options: e.g.,  -dump-ast-details -analyse-time-details
if [ -f "optitrust_flags.sh" ]; then
  source optitrust_flags.sh
fi


# NOTE: this could be removed if we use the VScode command to run "make optitrust" first.
# Run make update in working folder if requested
if [ "${RECOMPILE_OPTITRUST}" = "recompile_optitrust_yes" ]; then
  echo "recompile lib"
  make optitrust
  OUT=$?
  if [ ${OUT} -ne 0 ]; then
    echo "Could not compile lib"  >> /dev/stderr
    exit 1
  fi
fi


PROG="${FILEBASE}_with_lines.byte"

# First we create the source code for the transformation program
# ---DEPRECATED:
# ocaml ${VSCODE}/add_exit.ml -file "${FILEBASE}.ml" -line ${LINE}
#sed 's/show/myshow/' "${FILEBASE}.ml" > "${FILEBASE}_with_exit.byte"

# From "${FILEBASE}.ml", create ""{FILEBASE}_with_lines.ml" by inserting
# [~lines:__LINE__]   in the relevant places, and interpreting '!!' and '!!!'


# LATER: move this to a separate script
# Replace "!!^" with "Trace.check_exit_and_step ~line:__LINE__ ~reparse:true ()"
# Replace "!!!" with "Trace.check_exit_and_step ~line:__LINE__ ~ris_small_step:false ()"
# Replace "!!" with "Trace.check_exit_and_step ~line:__LINE__ ()"
sed 's/^\([[:space:]]*\)show /\1show ~line:__LINE__ /;s/\!\!\^/Trace.check_exit_and_step ~line:__LINE__ ~reparse:true ();/;s/!!!/Trace.check_exit_and_step ~line:__LINE__ ~is_small_step:false ();/;s/!!/Trace.check_exit_and_step ~line:__LINE__ ();/' "${FILEBASE}.ml" > "${FILEBASE}_with_lines.ml"
 # DEBUG: cat "${FILEBASE}_with_lines.ml"; exit 0

# LATER: add_exit should also introduce special commands for figuring out the line of the command that executes

# Second, we compile that transformation program
# DEPRECATED
# ocamlbuild -quiet -r -pkgs clangml,refl,pprint,str,optitrust "${FILEBASE}_with_exit.byte"
# TODO(Anton): replace this line with a dune command that uses directly /src/src files instead
# of the installed package; only consider ${FILEBASE}.ml from local folder


# TODO: check that PROG is also more recent than the optitrust library
if [[ "${FILEBASE}.ml" -nt "${PROG}" ]] || [[ "${FILEBASE}.cpp" -nt "${PROG}" ]]; then
  # echo FILE1 is newer than FILE2
  PROGNEEDSREBUILD="needsrebuild"
fi


if [ "${RECOMPILE_OPTITRUST}" = "recompile_optitrust_yes" ] || [ "${PROGNEEDSREBUILD}" = "needsrebuild" ]; then
  ocamlbuild -tag debug -quiet -r -pkgs clangml,refl,pprint,str,optitrust ${PROG}
fi

# LATER: capture the output error message
# so we can do the postprocessing on it


OUT=$?
if [ ${OUT} -ne 0 ];then
  echo "Could not compile file"  >> /dev/stderr
  exit 1
fi

# Third, we execute the transformation program, obtain "${FILEBASE}_before.cpp" and "${FILEBASE}_after.cpp
# Activate the backtrace
OCAMLRUNPARAM=b ./${PROG} -exit-line ${LINE} ${OPTIONS} ${OPTIONS2} ${FLAGS}

# DEBUG: echo "cd ${DIRNAME}; ./${PROG} -exit-line ${LINE} ${OPTIONS}"
# DEPREACTED | tee stdoutput.txt

OUT=$?
if [ ${OUT} -ne 0 ];then
  #echo "Error executing the script:"
  #echo "  cd ${DIRNAME}; ./${PROG} -exit-line ${LINE} ${OPTIONS} ${OPTIONS2}"
  exit 1
fi

# Fourth, we vizualize a result or a diff
# We need to cd to ${VSCODE} folder because that's how the scripts know the path to .vscode

cd ${VSCODE}


if [ "${VIEW}" = "view_diff" ] || [ "${VIEW}" = "view_diff_enc" ]; then

  if [ "${VIEW}" = "view_diff_enc" ]; then
    DIFFFOR="enc"
  fi

  ./open_diff.sh ${DIRNAME} ${FILEBASE} ${DIFFFOR} &

elif [ "${VIEW}" = "view_result" ]; then

  ./open_result.sh ${DIRNAME} ${FILEBASE} &

else

  echo "invalid VIEW argument"  >> /dev/stderr
  exit 1

fi

exit
