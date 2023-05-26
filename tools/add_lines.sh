#!/bin/sh

# Usage: add_lines.sh input.ml output.ml

# Replace "bigstep \"foo\"" with "Trace.open_bigstep ~line:__LINE__ \"foo\""
# Replace "!!" with "Trace.open_smallstep ~line:__LINE__ ()"
# Replace "!!!" with "Trace.open_smallstep ~line:__LINE__ ~reparse:true ()"


INPUT_FILE=$1
OUTPUT_FILE=$2

sed 's/^\([[:space:]]*\)show /\1show ~line:__LINE__ /
s/bigstep/Trace.open_bigstep ~line:__LINE__ /
s/!!!/Trace.open_smallstep ~line:__LINE__ ~reparse:true ();/
s/!!/Trace.open_smallstep ~line:__LINE__ ();/' ${INPUT_FILE} > ${OUTPUT_FILE}
