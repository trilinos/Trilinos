#!/bin/bash -e

#  Build a generic version of <Project>BuildQuickRef for general
#  documentation.

source ../utils/gen_doc_utils.sh

generate_git_version_file

echo
echo "Generating HTML and PDF files ..."
echo

make
