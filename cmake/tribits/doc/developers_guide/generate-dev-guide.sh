#!/bin/bash -e

# This script is used to generate the the TribitsDevelopersGuide.(html,pdf)
# files using a new script in TriBITS for that purpose.  You just run it from
# this directory as:
#
#   cd <this-dir>
#   ./generate-dev-guide.sh
#
# If you want to override what utils are used to generate the files, you can
# pass in, for exmaple:
#
#   ./generate-dev-guide.sh \
#      --generate-html=rst2html.py --generate-latex=rst2latex.py
#
# Note that if you are debugging the parsing, you can disable the generation
# of the latex and pdf by setting:
#
#    --generate-latex=
#
# This script also automatically extracts detailed TriBITS documentation from
# the *.cmake files using extract_rst_cmake_doc.py (which works kind of like
# doxygen).
#
# To see output from extract_rst_cmake_doc.py just run the script as:
#
#   env EXTRACT_RST_CMAKE_DOC_EXTRA_ARGS=--do-trace \
#      ./generate-dev-guide.sh [other args]
#
# NOTE: If you see rst2html or rst2latex errors for the file
# TribitsDeveloeprsGuilde.rst with line numbers that don't seem to make sense,
# this is likley due to the include of TribitsDetailedMacroFunctionDoc.rst.
# To adjust the line numbers, subtract the line number of the include for
# TribitsDetailedMacroFunctionDoc.rst in TribitsDevelopersGuilde.rst form the
# line number given in the output and that will be the line number in the
# TribitsDevelopersGuilde.rst file.  You can then match that up with the
# original text in the *.cmake file that this came from for the given macro or
# function.
#
# Enjoy!

ARGS=$@

echo
echo "Extracting TriBITS documentation from *.cmake files ..."
echo
../../python/extract_rst_cmake_doc.py \
  --extract-from=../../package_arch/,../../utils/ \
  --rst-file-pairs=TribitsDetailedMacroFunctionDocTemplate.rst:TribitsDetailedMacroFunctionDoc.rst \
  $EXTRACT_RST_CMAKE_DOC_EXTRA_ARGS



echo
echo "Generating HTML and PDF files ..."
echo
../../python/generate-docutils-output.py \
  --file-base=TribitsDevelopersGuide \
  $ARGS

# NOTE: This above invocation by default overrides the DocUtils html and latex
# generators to use the versions without the *.py prefix.  When you install
# docutils from source, you get the *.py prefix.  When you install using the
# yum package on a Linux machine, you don't get the *.py prefix.  What a pain!
