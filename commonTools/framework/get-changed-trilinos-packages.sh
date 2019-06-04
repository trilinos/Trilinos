#!/bin/bash -e
#
# Usage:
#
#   get-changed-trilinos-packages.sh <git-commit-from> <git-commit-to> \
#     <package-enables-cmake-out>
#
# This script takes a range of git commits <git-commit-from>..<git-commit-to>
# and then generates a CMake fragment file <package-enables-cmake-out> which
# provides the set of enables of Trilinos packages needed to test the changed
# files.
#
# For example, to generate a file for the set of enables to test changes in
# the current version of Trilinos w.r.t. to develop branch, one would do:
#
#   $ cd <build-dir>/
#
#   $ $TRILINOS_DIR/commonTools/test/utilities/get-changed-trilinos-packages.sh \
#     origin/develop HEAD packageEnables.cmake
#
#   $ cmake -C packageEnables.cmake [other cache vars] $TRILINOS_DIR
#
# By default, the Trilinos git repo base directory is determined automatically
# when this script is run directly from it's home location in Trilinos.
# Otherwise, the env var TRILINOS_DIR must be set to point to the base
# Trilinos git repo before calling this script.
#
# As a bi-product of running this script it will create the following files in
# the current working directory:
#
# * TribitsDumpDepsXmlScript.log: Log file from running TribitsDumpDepsXmlScript.cmake
# * TrilinosPackageDependencies.xml: Generated XML package list and dependencies
# * changed-files.txt: List of files changed between git repo versions
#
# The file <package-enables-cmake-out> can be created in any directory by
# giving the relative or absolute path.

#
# Get command-line arguments
#

GIT_COMMIT_FROM=$1
GIT_COMMIT_TO=$2
CMAKE_PACKAGE_ENABLES_OUT=$3

if [ "$GIT_COMMIT_FROM" == "" ] ; then
  echo "ERROR: Must specify first argument <git-commit-from>!"
  exit 1
fi

if [ "$GIT_COMMIT_TO" == "" ] ; then
  echo "ERROR: Must specify second argument <git-commit-to>!"
  exit 2
fi

if [ "$CMAKE_PACKAGE_ENABLES_OUT" == "" ] ; then
  echo "ERROR: Must specify third argument <package-enables-cmake-out>!"
  exit 1
fi

echo
echo "***"
echo "*** Generating set of Trilinos enables given modified packages from"
echo "*** git commit ${GIT_COMMIT_FROM} to ${GIT_COMMIT_TO}"
echo "***"
echo

#
# Determine TRILINOS_DIR
#

if [ "$TRILINOS_DIR" == "" ] ; then
  # Grab from the symlink (only works on Linux)
  _ABS_FILE_PATH=`readlink -f $0` || \
   echo "Could not follow symlink to set TRILINOS_DIR!"
  if [ "$_ABS_FILE_PATH" != "" ] ; then
    _SCRIPT_DIR=`dirname $_ABS_FILE_PATH`
    TRILINOS_DIR=$_SCRIPT_DIR/../..
  fi
fi

if [ "$TRILINOS_DIR" == "" ] ; then
  echo "ERROR: Cannot determine TRILINOS_DIR!  Please set env var TRILINOS_DIR!"
  exit 4
fi

echo "TRILINOS_DIR=$TRILINOS_DIR"

# Allow a different source tree for the Trilinos scripts
if [ "$TRILINOS_SCRIPTS_DIR" == "" ] ; then
  TRILINOS_SCRIPTS_DIR=${TRILINOS_DIR}
fi

echo "TRILINOS_SCRIPTS_DIR=$TRILINOS_SCRIPTS_DIR"

ORIG_CWD=$PWD

# Allow override of TriBITS for testing purposes
if [ "${GCTP_TRIBITS_DIR_OVERRIDE}" != "" ] ; then
  TRIBITS_DIR=$GCTP_TRIBITS_DIR_OVERRIDE
else
  TRIBITS_DIR=$TRILINOS_SCRIPTS_DIR/cmake/tribits
fi
echo "TRIBITS_DIR=$TRIBITS_DIR"

echo
echo "A) Generate the Trilinos Packages definition and depencencies XML file"
echo

cmake \
  -D Trilinos_DEPS_XML_OUTPUT_FILE=TrilinosPackageDependencies.xml \
  -P $TRIBITS_DIR/ci_support/TribitsDumpDepsXmlScript.cmake \
  &> TribitsDumpDepsXmlScript.log

echo "Wrote the file 'TrilinosPackageDependencies.xml'"

echo
echo "B) Get the set of changed files"
echo

# Allow override of git for unit testing purposes
if [ "${GCTP_GIT_OVERRIDE}" != "" ] ; then
  GIT_EXEC=${GCTP_GIT_OVERRIDE}
else
  GIT_EXEC=git
fi

CHANGED_FILES_FILE=$ORIG_CWD/changed-files.txt

cd $TRILINOS_DIR/
echo "Current directory: $PWD"
echo
echo "$GIT_EXEC diff --name-only ${GIT_COMMIT_FROM}..${GIT_COMMIT_TO} > ${CHANGED_FILES_FILE}"
$GIT_EXEC diff --name-only ${GIT_COMMIT_FROM}..${GIT_COMMIT_TO} > ${CHANGED_FILES_FILE}
cd $ORIG_CWD/
echo
echo "Wrote file 'changed-files.txt'"
echo
echo "Current directory: $PWD"

echo
echo "C) Get the unfiltered list of changed Trilinos packages (including 'ALL_PACKAGES')"
echo
CHANGED_PACKAGES_FULL_LIST=`$TRIBITS_DIR/ci_support/get-tribits-packages-from-files-list.py \
  --deps-xml-file=TrilinosPackageDependencies.xml \
  --files-list-file=changed-files.txt`
echo "CHANGED_PACKAGES_FULL_LIST='$CHANGED_PACKAGES_FULL_LIST'"

echo
echo "D) Filter list of changed packages to get only the PT packages"
echo
CHANGED_PACKAGES_PT_LIST=`$TRIBITS_DIR/ci_support/filter-packages-list.py \
  --deps-xml-file=TrilinosPackageDependencies.xml \
  --input-packages-list=$CHANGED_PACKAGES_FULL_LIST \
  --keep-test-test-categories=PT`
echo "CHANGED_PACKAGES_PT_LIST='$CHANGED_PACKAGES_PT_LIST'"

echo
echo "E) Generate the *.cmake enables file"
echo

echo "
MACRO(PR_ENABLE_BOOL  VAR_NAME  VAR_VAL)
  MESSAGE(\"-- Setting \${VAR_NAME} = \${VAR_VAL}\")
  SET(\${VAR_NAME} \${VAR_VAL} CACHE BOOL \"Set in $CMAKE_PACKAGE_ENABLES_OUT\")
ENDMACRO()
" >  $CMAKE_PACKAGE_ENABLES_OUT

if [ "$CHANGED_PACKAGES_PT_LIST" != "" ] ; then
  echo "$CHANGED_PACKAGES_PT_LIST" | sed -n 1'p' | tr ',' '\n' | while read PKG_NAME ; do
    #echo $PKG_NAME
    echo "PR_ENABLE_BOOL(Trilinos_ENABLE_${PKG_NAME} ON)" >> $CMAKE_PACKAGE_ENABLES_OUT
  done
else
  echo "MESSAGE(\"-- NOTE: No packages are being enabled!\")" >> $CMAKE_PACKAGE_ENABLES_OUT
fi

echo "Wrote file '$CMAKE_PACKAGE_ENABLES_OUT'"
