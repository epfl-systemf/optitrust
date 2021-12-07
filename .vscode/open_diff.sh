#!/bin/bash

# How to test this script:
# cd src/.vscode
# assume test_split_before.cpp and test_split_after.cpp exist
# ./open_diff.sh ../unit_tests test_split

VSCODE=`pwd`
TOOLS_FOLDER="${VSCODE}/../tools"


DIRNAME=$1
FILEBASE=$2

# Work in the file directory
cd ${DIRNAME}

TARGET="${FILEBASE}_diff.html"

#echo "Generating ${TARGET}"

# Compute title
TITLESTR="OptiTrust_Diff_${FILEBASE}"


# Read options from optitrust_flags.sh
# options: e.g., FLAGS="-dump-ast-details -analyse-time-details"
if [ -f "optitrust_flags.sh" ]; then
  source optitrust_flags.sh
fi
if [ "${DIFF_ENC}" = "true" ]; then
  ENCOPT="_enc"
fi

# Compute diff
DIFFCODE=`git diff --ignore-all-space --no-index -U10 ${FILEBASE}_before${ENCOPT}.cpp ${FILEBASE}_after${ENCOPT}.cpp | base64 -w 0`
DIFFSTR="var diffString = window.atob(\"${DIFFCODE}\");"

# Take templace and substitute ${TOOLS_FOLDER}, ${INSERT_TITLE}, and ${INSERT_DIFF}
TEMPLATE="${TOOLS_FOLDER}/diff_template.html"

cp ${TEMPLATE} ${TARGET}

sed -i "s#{INSERT_TITLE}#${TITLESTR}#g;s#{TOOLS_FOLDER}#${TOOLS_FOLDER}#g;s#{INSERT_DIFF}#${DIFFSTR}#g" ${TARGET}
# Note: there seems to be an issue if performing the sed one after the other...
#echo "Generated ${TARGET}"

# Open the browser with the target file
cd ${VSCODE}
./open_in_browser.sh ${DIRNAME}/${TARGET} "${TITLESTR}"
