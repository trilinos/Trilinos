#!/bin/bash

#
# Symlink into some scratch directory and then run as:
#
#   ./ctest-s-local-test-driver.sh [build-name-1 build-name-2 ...]
#
# Pass in '--help' or '-h' to print help.
#

CTEST_S_LOCAL_DRIVER_HELP_STR="Drive builds locally and submit to CDash

Usage: ./ctest-s-local-test-driver.sh <build-name-keys0> <build-name-keys1> ...

To use, symlink into some scratch directory and then run as:

  ./ctest-s-local-test-driver.sh <build-name-1> <build-name-2> ...

If no build names are not given, then the help message is printed.  If 'all'
is given, then all of the builds supported for the current sytem are run as
specified in the Trilinos/cmake/std/atdm/<system_name>/all_supported_builds.sh
file.

The buld names <build-name-keysi> by default must match the the name of the
builds listed under:

  Trilinos/cmake/ctest/drivers/atdm/<system_name>/drivers/

By default, <build-name-keysi> is the partial build name after the build-name
prefix to form the full build name:

  ${ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX}<build-name-keysi>

(where ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX is defined in the
Trilinos/cmake/std/atdm/<system_name>/all_supported_builds.sh file) and the
full driver script name:

  Trilinos/cmake/ctest/drivers/atdm/<system_name>/drivers/
    ${ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX}<build-name-keysi>.sh

For example, ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX=Trilinos-atdm-waterman- and
<build-name-keysi>=cuda-9.2-release-debug gives the full build name
'Trilinos-atdm-waterman-cuda-9.2-release-debug'.

However, if env var:

  ATDM_CTEST_S_USE_FULL_BUILD_NAME=1

is set, then <build-name-keysi> must match the full build name
(e.g. <build-name-keysi>=Trilinos-atdm-waterman_cuda-9.2_fpic_opt) giving the
full script name:

  Trilinos/cmake/ctest/drivers/atdm/<system_name>/drivers/
    <build-name-keysi>.sh

Tail -f the generated files <full_build_name>/smart-jenkins-driver.out to see
details of each run (e.g. <full_build_name> =
Trilinos-atdm-<system_name>-gnu-openmp-opt).

To select the default env to load instead of 'default', use:

  env ATDM_CHT_DEFAULT_ENV=<system_name>-default \\
  ./ctest-s-local-test-driver.sh <build-name-1> >build-name-2> ...

(For example, this is needed for the 'cee-rhel6' system to set
ATDM_CHT_DEFAULT_ENV=cee-rhel6-default, otherwise the 'sems-rhel6' env will be
selected.)

To control the list of packages tests, not rebuild from scratch, and not
submit, use (for example):

  env \\
    Trilinos_PACKAGES=<pkg0>,<pkg1>,... \\
    CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE \\
    CTEST_DO_SUBMIT=OFF \\
  ./ctest-s-local-test-driver.sh <build-name-1> >build-name-2> ...

Other options that are good to set sometimes include:

  CTEST_DO_CONFIGURE=OFF
  CTEST_DO_BUILD=OFF
  CTEST_BUILD_FLAGS='<flag0> <flag1> ...'
  CTEST_DO_TEST=OFF
  CTEST_PARALLEL_LEVEL=<N>
  CTEST_DO_SUBMIT=OFF
"

if [[ "$@" == "" ]] || [[ "$@" == "-h" ]] ||  [[ "$@" == "--help" ]]; then
  echo "$CTEST_S_LOCAL_DRIVER_HELP_STR"
  exit 0
fi

#
# Sound off
#

echo
echo "***"
echo "*** $0"
echo "***"
echo

#
# Get the location of the base Trilinos directory
#

if [ "$ATDM_TRILINOS_DIR" == "" ] ; then
  # Grab from the symlink (only works on Linux)
  _ABS_FILE_PATH=`readlink -f $0` || \
   echo "Could not follow symlink to set TRILINOS_DIR!"
  if [ "$_ABS_FILE_PATH" != "" ] ; then
    export STD_ATDM_DIR=`dirname $_ABS_FILE_PATH`
    export ATDM_TRILINOS_DIR=`readlink -f $STD_ATDM_DIR/../../..`
  fi
fi

export STD_ATDM_DIR="${ATDM_TRILINOS_DIR}/cmake/std/atdm"

echo "ATDM_TRILINOS_DIR = '$ATDM_TRILINOS_DIR'"

if [ "$ATDM_TRILINOS_DIR" == "" ] ; then
  echo "ERROR: Cannot determine TRILINOS_DIR (you must be on a non-Linux system or you must have copied the script instead of symlinking it as per instructions)."
  exit 1
fi

#
# Load a default env for the system
#

if [ "$ATDM_CHT_DEFAULT_ENV" == "" ] ; then
  ATDM_CHT_DEFAULT_ENV=default
fi
#echo "ATDM_CHT_DEFAULT_ENV = ${ATDM_CHT_DEFAULT_ENV}"

echo
echo "Load some env to get python, cmake, etc ..."
echo
source $STD_ATDM_DIR/load-env.sh ${ATDM_CHT_DEFAULT_ENV}
# NOTE: Above, it does not matter which env you load.  Any of them will
# provide the right python, cmake, etc.

#
# Get the list of builds to run
# 

# Must get ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX
source $STD_ATDM_DIR/$ATDM_CONFIG_SYSTEM_NAME/all_supported_builds.sh
#echo "ATDM_CONFIG_ALL_SUPPORTED_BUILDS = '${ATDM_CONFIG_ALL_SUPPORTED_BUILDS[@]}'"

ATDM_ARRAY_OF_BUILDS=$@
if [ "${ATDM_ARRAY_OF_BUILDS}" == "all" ] ; then
  ATDM_ARRAY_OF_BUILDS=${ATDM_CONFIG_ALL_SUPPORTED_BUILDS[@]}
fi

echo
echo "Running builds:"
for build_name in ${ATDM_ARRAY_OF_BUILDS[@]} ; do
  echo "    ${build_name}"
done

#
# Run the builds using the ctest -S driver script
#

BASEDIR=$PWD

ln -sf ${ATDM_TRILINOS_DIR} .

for build_name_body in ${ATDM_ARRAY_OF_BUILDS[@]} ; do

  if [ "${ATDM_CTEST_S_USE_FULL_BUILD_NAME}" == "1" ] ; then
    build_name="${build_name_body}"
  else
    build_name="${ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX}${build_name_body}"
  fi

  echo
  echo "Running Jenkins driver ${build_name}.sh ..."

  # Set up the directory for this build case

  if [ ! -e ${build_name} ] ; then
    echo
    echo "    Creating directory: ${build_name}"
    mkdir ${build_name}
  fi
  cd ${BASEDIR}/${build_name}
  #pwd

  #echo "Creating symlink: Trilinos" 
  ln -sf ${ATDM_TRILINOS_DIR} .

  # Set up the SRC_AND_BUILD dir for Trilinos already cloned

  if [ ! -e SRC_AND_BUILD ] ; then
    echo
    echo "    Creating directory: SRC_AND_BUILD" 
    mkdir SRC_AND_BUILD
  fi

  #echo "Creating symlink: SRC_AND_BUILD/Trilinos" 
  cd SRC_AND_BUILD/
  ln -sf ${ATDM_TRILINOS_DIR} .
  cd ..

  echo
  echo "    See log file ${build_name}/smart-jenkins-driver.out"

  time env \
    JOB_NAME=${build_name} \
    WORKSPACE=$PWD \
    CTEST_TEST_TYPE=Experimental \
    CTEST_DO_UPDATES=OFF \
  ${ATDM_TRILINOS_DIR}/cmake/ctest/drivers/atdm/smart-jenkins-driver.sh \
    &> smart-jenkins-driver.out

  cd ${BASEDIR}

done
