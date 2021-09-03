#!/bin/bash -e
#
# Build all of the TriBITS-related documentation
#
# To build this documentation, from any directory, run:
#
#   <base-dir>/build_docs.sh
#
# To only produce all of the input *.rst files but skip the final generation
# of the HTML and other output files, run:
#
#   <base-dir>/build_docs.sh --skip-final-generation
#

#
# A) Parse command-line arguments
#

skip_final_generation=0
skip_final_generation_arg=

while (( "$#" )); do
  case "$1" in
    --skip-final-generation)
      skip_final_generation=1
      skip_final_generation_arg=--skip-final-generation
      shift
      ;;
    *)
      echo "Error: The argument '$1' is not supported!"
      exit 1
      ;;
  esac
done

#
# B) CD into the tribits/doc directory
#

if [ "$TRIBITS_BASE_DIR" == "" ] ; then
  _ABS_FILE_PATH=`readlink -f $0`
  _BASE_DIR=`dirname $_ABS_FILE_PATH`
fi

cd $_BASE_DIR

#
# C) Generate the documentation
#

echo
echo "***"
echo "*** Generating Users and Maintainers guides ..."
echo "***"
echo

cd $_BASE_DIR/guides/
./generate-guide.sh all ${skip_final_generation_arg}
cd -
if [[ "${skip_final_generation}" == "0" ]] ; then
  echo
  echo "See generated files:"
  echo
  ls $_BASE_DIR/guides/users_guide/TribitsUsersGuide.*
  ls $_BASE_DIR/guides/maintainers_guide/TribitsMaintainersGuide.*
fi

echo
echo "***"
echo "*** Generating TribitsBuildReference.[rst,html,pdf] ..."
echo "***"
echo

cd $_BASE_DIR/build_ref
./create-build-ref.sh ${skip_final_generation_arg}
if [[ "${skip_final_generation}" == "0" ]] ; then
  echo
  echo "See generated files:"
  echo
  ls $_BASE_DIR/build_ref/TribitsBuildReference.*
fi

#echo
#echo "***"
#echo "*** Generating TribitsOverview.pdf ..."
#echo "***"
#echo
#
#cd $_BASE_DIR/overview
#source source_set_env 
#make
#echo
#echo "See generated files:"
#echo
#ls $_BASE_DIR/overview/TribitsOverview.pdf
#
#echo
#echo "***"
#echo "*** Generating TribitsLifecycleModel.pdf ..."
#echo "***"
#echo
#
#cd $_BASE_DIR/lifecycle_model
#make
#echo
#echo "See generated files:"
#echo
#ls $_BASE_DIR/lifecycle_model/TribitsLifecycleModel.pdf
