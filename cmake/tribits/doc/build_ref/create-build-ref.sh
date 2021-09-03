#!/bin/bash -e
#
# Build a generic version of <Project>BuildReference for general
# documentation.
#
# Run as:
#
#   cd <this_dir>/
#   ./create-build-ref.sh
#
# You can skip the final generation of the output HTML and other files by
# calling:
#
#   ./create-build-ref.sh --skip-final-generation
# 

skip_final_generation=0

while (( "$#" )); do
  case "$1" in
    --skip-final-generation)
      skip_final_generation=1
      shift
      ;;
    *)
      echo "Error: The argument '$1' is not supported!"
      exit 1
      ;;
  esac
done

source ../utils/gen_doc_utils.sh

generate_git_version_file

if [[ "${skip_final_generation}" == "0" ]] ; then
  echo
  echo "Generating HTML and PDF files ..."
  echo
  make
else
  echo
  echo "Only generating the file TribitsBuildReference.rst on request!"
  echo
  make TribitsBuildReference.rst
fi
