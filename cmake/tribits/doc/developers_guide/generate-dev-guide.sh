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
#   $ env TRIBITS_DEV_GUIDE_EXTRACT_RST_CMAKE_DOC_EXTRA_ARGS=--do-trace \
#       ./generate-dev-guide.sh [other args]
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
# NOTE: To skip the extraction of the documentation from the *.cmake files,
# just sent the env TRIBITS_DEV_GUIDE_varaible SKIP_DOCUMENTATION_EXTRACTION
# as:
#
#   $ env TRIBITS_DEV_GUIDE_SKIP_DOCUMENTATION_EXTRACTION=1 \
#      ./generate-dev-guilde.sh
#
# That will result in the generated files TribitsMacroFunctionDoc.rst and
# UtilsMacroFunctionDoc.rst being left as is.  This would be useful to speed
# up builds (but it is very fast) but is more useful when spell checking and
# editing the documentation.  This speeds up the editing process and then the
# updated documentation can be copied back into the *.cmake files of origin.
#
# Enjoy!

ARGS=$@

source source_set_env


if [ "$TRIBITS_DEV_GUIDE_SKIP_DOCUMENTATION_EXTRACTION" == "" ] ; then

  echo
  echo "Extracting TriBITS documentation from *.cmake files ..."
  echo
  ../../python_utils/extract_rst_cmake_doc.py \
    --extract-from=../../core/package_arch/,../../core/utils/,../../ctest_driver/ \
    --rst-file-pairs=TribitsMacroFunctionDocTemplate.rst:TribitsMacroFunctionDoc.rst,UtilsMacroFunctionDocTemplate.rst:UtilsMacroFunctionDoc.rst \
    $TRIBITS_DEV_GUIDE_EXTRACT_RST_CMAKE_DOC_EXTRA_ARGS
  
fi

if [ "$TRIBITS_DEV_GUIDE_SKIP_OTHER_EXTRACTION" == "" ] ; then

  echo
  echo "Generating list of Standard TriBITS TPLs ..."
  echo
  ls -w 1 ../../core/std_tpls/ &> TribitsStandardTPLsList.txt

  echo
  echo "Generating Directory structure of TribitsHelloWorld ..."
  echo
  ../../python_utils/tree.py -f -c -x ../../examples/TribitsHelloWorld/ \
    &> TribitsHelloWorldDirAndFiles.txt

  echo
  echo "Generating output for 'checkin-test.py --help' ..."
  echo
  ../../ci_support/checkin-test.py --help &> checkin-test-help.txt

  echo
  echo "Generating output for 'gitdist --help' ..."
  echo
  ../../python_utils/gitdist --help &> gitdist-help.txt

  echo
  echo "Generating output for 'snapshot-dir.py --help' ..."
  echo
  env SNAPSHOT_DIR_DUMMY_DEFAULTS=1 ../../python_utils/snapshot-dir.py --help \
   &> snapshot-dir-help.txt

fi


echo
echo "Generating HTML and PDF files ..."
echo
../../python_utils/generate-docutils-output.py \
  --file-base=TribitsDevelopersGuide \
  --generate-latex-options="--stylesheet-path=rst2latex.tex" \
  $ARGS
