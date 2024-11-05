#!/bin/bash

# This script is used to generate the TribitsUsersGuide.(html,pdf) and
# TribitsMaintainersGuide.html files.  You just run it from this directory as:
#
#   cd <this-dir>
#   ./generate-guide.sh all
#
# To just generate the users guide, run:
#
#   ./generate-guide.sh ug
#
# To just generate the maintainers guide, run:
#
#   ./generate-guide.sh mg
#
# This script automatically extracts detailed TriBITS documentation from the
# *.cmake files using the tool extract_rst_cmake_doc.py (which works kind of
# like doxygen).  To see output from extract_rst_cmake_doc.py just run the
# script as:
#
#   $ env TRIBITS_DEV_GUIDE_EXTRACT_RST_CMAKE_DOC_EXTRA_ARGS=--do-trace \
#       ./generate-guide.sh all
#
# NOTE: If you see rst2html or rst2latex errors in the files
# TribitsUsersGuide.rst or TribitsMaintainersGuide.rst with line numbers that
# don't seem to make sense, this is likely due to the include of
# TribitsDetailedMacroFunctionDoc.rst.  To adjust the line numbers, subtract
# the line number of the include for TribitsDetailedMacroFunctionDoc.rst in
# TribitsUsersGuide.rst or TribitsMaintainersGuide.rst from the line number
# given in the output and that will be the line number in the file
# TribitsUsersGuide.rst or TribitsMaintainersGuide.rst.  You can then match
# that up with the original text in the *.cmake file that this came from for
# the given macro or function.
#
# NOTE: To skip the extraction of the documentation from the *.cmake files,
# just sent the env TRIBITS_DEV_GUIDE_SKIP_DOCUMENTATION_EXTRACTION variable
# as:
#
#   $ env TRIBITS_DEV_GUIDE_SKIP_DOCUMENTATION_EXTRACTION=1 \
#      ./generate-guide.sh
#
# That will result in the generated files TribitsMacroFunctionDoc.rst and
# UtilsMacroFunctionDoc.rst being left as is.  This would be useful to speed
# up builds (it is very fast) but is more useful when spell checking and
# editing the documentation.  This speeds up the editing process and then the
# updated documentation can be copied back into the *.cmake files of origin.
#
# Enjoy!
#


#
# Parse command-line arguments
#

generate_maintainers_guide=0
generate_users_guide=0
skip_final_generation=0

while (( "$#" )); do
  case "$1" in
    mg|maintainers_guide)
      generate_maintainers_guide=1
      shift
      ;;
    ug|users_guide)
      generate_users_guide=1
      shift
      ;;
    all)
      generate_maintainers_guide=1
      generate_users_guide=1
      shift
      ;;
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


#
# Functions
#

source ../utils/gen_doc_utils.sh


function generate_gitdist_dist_help_topic {
  help_topic_name=$1
  ../../python_utils/gitdist --dist-help=$help_topic_name &> gitdist-dist-help-$help_topic_name.txt.tmp
  update_if_different  gitdist-dist-help-$help_topic_name.txt  tmp
}


function tribits_extract_rst_cmake_doc {
  dir=$1
  extra_args=$2

  cd ${dir}
  echo $PWD

  if [ "$TRIBITS_DEV_GUIDE_SKIP_DOCUMENTATION_EXTRACTION" == "" ] ; then

    echo
    echo "Extracting TriBITS documentation from *.cmake files ..."
    echo
    ../../../python_utils/extract_rst_cmake_doc.py \
      --extract-from=../../../ctest_driver/,../../../ci_support/,../../../core/package_arch/,../../../core/test_support/,../../../core/utils/ \
      --rst-file-pairs=../TribitsMacroFunctionDocTemplate.rst:TribitsMacroFunctionDoc.rst.tmp,../UtilsMacroFunctionDocTemplate.rst:UtilsMacroFunctionDoc.rst.tmp,../TribitsSystemMacroFunctionDocTemplate.rst:TribitsSystemMacroFunctionDoc.rst.tmp \
      ${extra_args} \
      --file-name-path-base-dir=../../.. \
      $TRIBITS_DEV_GUIDE_EXTRACT_RST_CMAKE_DOC_EXTRA_ARGS

    update_if_different  TribitsMacroFunctionDoc.rst  tmp
    update_if_different  UtilsMacroFunctionDoc.rst  tmp
    update_if_different  TribitsSystemMacroFunctionDoc.rst  tmp

  fi

  cd -

}


function tribits_extract_other_doc {

  if [ "$TRIBITS_DEV_GUIDE_SKIP_OTHER_EXTRACTION" == "" ] ; then

    echo
    echo "Generating list of Standard TriBITS TPLs ..."
    echo
    ls -w 1 ../../core/std_tpls/ &> TribitsStandardTPLsList.txt.tmp
    update_if_different  TribitsStandardTPLsList.txt  tmp

    echo
    echo "Generating list of Common TriBITS TPLs ..."
    echo
    ls -w 1 ../../common_tpls/ &> TribitsCommonTPLsList.txt.tmp
    update_if_different  TribitsCommonTPLsList.txt  tmp

    echo
    echo "Generating Directory structure of TribitsHelloWorld ..."
    echo
    ../../python_utils/tree.py -f -c -x ../../examples/TribitsHelloWorld/ \
      &> TribitsHelloWorldDirAndFiles.txt.tmp
    update_if_different  TribitsHelloWorldDirAndFiles.txt  tmp

    echo
    echo "Generating TribitsExampleProject2/Package1 CMakeList file variants ..."
    echo

    cat ../../examples/TribitsExampleProject2/packages/package1/CMakeLists.raw.cmake \
      | grep -v EnableTribitsTestSupport \
      | grep -v GeneratePackageConfigFileForBuildDir \
      &> TribitsExampleProject2_Package1_CMakeLists.raw.external.cmake.tmp
    update_if_different TribitsExampleProject2_Package1_CMakeLists.raw.external.cmake tmp

    cat ../../examples/TribitsExampleProject2/packages/package1/CMakeLists.raw.cmake \
      | grep -v EnableTribitsTestSupport \
      &> TribitsExampleProject2_Package1_CMakeLists.raw.internal.cmake.tmp
    update_if_different TribitsExampleProject2_Package1_CMakeLists.raw.internal.cmake tmp

    cat ../../examples/TribitsExampleProject2/packages/package1/CMakeLists.raw.cmake \
      | grep -v "that TriBITS does automatically" \
      | grep -v DefineAllLibsTarget \
      | grep -v GeneratePackageConfigFileForBuildDir \
      | grep -v GeneratePackageConfigFileForInstallDir \
      &> TribitsExampleProject2_Package1_CMakeLists.raw.tribits_test.cmake.tmp
    update_if_different \
      TribitsExampleProject2_Package1_CMakeLists.raw.tribits_test.cmake tmp

    echo
    echo "Generating output for 'checkin-test.py --help' ..."
    echo
    ../../ci_support/checkin-test.py --help &> checkin-test-help.txt.tmp
    update_if_different  checkin-test-help.txt  tmp

    echo
    echo "Generating output for 'gitdist --help' and '--dist-help=<topic>' ..."
    echo
    ../../python_utils/gitdist --help &> gitdist-help.txt.tmp
    update_if_different  gitdist-help.txt  tmp
    generate_gitdist_dist_help_topic overview
    generate_gitdist_dist_help_topic repo-selection-and-setup
    generate_gitdist_dist_help_topic dist-repo-status
    generate_gitdist_dist_help_topic repo-versions
    generate_gitdist_dist_help_topic aliases
    generate_gitdist_dist_help_topic default-branch
    generate_gitdist_dist_help_topic move-to-base-dir
    generate_gitdist_dist_help_topic usage-tips
    generate_gitdist_dist_help_topic script-dependencies
    generate_gitdist_dist_help_topic all

    echo
    echo "Generating output for 'clone_extra_repos.py --help' ..."
    echo
    ../../ci_support/clone_extra_repos.py --help \
      &> clone_extra_repos-help.txt.tmp
    update_if_different  clone_extra_repos-help.txt  tmp

    echo
    echo "Generating output for 'snapshot-dir.py --help' ..."
    echo
    env SNAPSHOT_DIR_DUMMY_DEFAULTS=1 ../../python_utils/snapshot-dir.py --help \
     &> snapshot-dir-help.txt.tmp
    update_if_different  snapshot-dir-help.txt  tmp

    echo
    echo "Generating output for 'is_checkin_tested_commit.py --help' ..."
    echo
    ../../ci_support/is_checkin_tested_commit.py --help \
      &> is_checkin_tested_commit.txt.tmp
    update_if_different  is_checkin_tested_commit.txt  tmp

    echo
    echo "Generating output for 'get-tribits-packages-from-files-list.py --help' ..."
    echo
    ../../ci_support/get-tribits-packages-from-files-list.py --help \
      &> get-tribits-packages-from-files-list.txt.tmp
    update_if_different  get-tribits-packages-from-files-list.txt  tmp

    echo
    echo "Generating output for 'get-tribits-packages-from-last-tests-failed.py --help' ..."
    echo
    ../../ci_support/get-tribits-packages-from-last-tests-failed.py --help \
      &> get-tribits-packages-from-last-tests-failed.txt.tmp
    update_if_different  get-tribits-packages-from-last-tests-failed.txt  tmp

    echo
    echo "Generating output for 'filter-packages-list.py --help' ..."
    echo
    ../../ci_support/filter-packages-list.py --help \
      &> filter-packages-list.txt.tmp
    update_if_different  filter-packages-list.txt  tmp

    echo
    echo "Generating output for 'install_devtools.py --help' ..."
    echo
    ../../devtools_install/install_devtools.py --help \
      &> install_devtools-help.txt.tmp
    update_if_different  install_devtools-help.txt  tmp

  fi

  if [ -e "../../../README.DIRECTORY_CONTENTS.rst" ] ; then
    echo
    echo "Copy TriBITS/README.DIRECTORY_CONTENTS.rst to TriBITS.README.DIRECTORY_CONTENTS.rst ..."
    echo
    cp ../../../README.DIRECTORY_CONTENTS.rst TriBITS.README.DIRECTORY_CONTENTS.rst.tmp
  else
    echo
    echo "TriBITS/README.DIRECTORY_CONTENTS.rst does not exist to copy!"
    echo
    touch TriBITS.README.DIRECTORY_CONTENTS.rst.tmp
  fi
  update_if_different  TriBITS.README.DIRECTORY_CONTENTS.rst  tmp

}


function make_final_doc_in_subdir {
  dir_name=$1
  if [[ "${skip_final_generation}" == "0" ]] ; then
    cd $dir_name
    echo $PWD
    time make
    cd - > /dev/null
  else
    echo
    echo "Skipping final generation of '${dir_name}' on request!"
    echo
  fi
}


#
# Executable code
#

generate_git_version_file

tribits_extract_other_doc

if [[ "${generate_users_guide}" == "1" ]] ; then

  echo
  echo "Generating HTML and PDF files for users guide ..."
  echo

  tribits_extract_rst_cmake_doc  users_guide

  make_final_doc_in_subdir  users_guide

fi

if [[ "${generate_maintainers_guide}" == "1" ]] ; then

  echo
  echo "Generating HTML file for maintainers_guide ..."
  echo

  tribits_extract_rst_cmake_doc  maintainers_guide --show-file-name-line-num

  make_final_doc_in_subdir  maintainers_guide

fi
