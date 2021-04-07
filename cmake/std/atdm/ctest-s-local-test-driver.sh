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

To run all of the supported builds, run with 'all':

  ./ctest-s-local-test-driver.sh all

which runs all of the supported builds listed in the file
Trilinos/cmake/std/atdm/<system_name>/all_supported_builds.sh.

If no commandline arguments are given, then then the list of supported builds
that can be selected from is printed.

If specifying the individual names <build-name-keysi> then the much match the
names of the driver scripts listed under:

  Trilinos/cmake/ctest/drivers/atdm/<system_name>/drivers/

By default, <build-name-keysi> is the partial build name after the build-name
prefix to form the full build name:

  ${ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX}<build-name-keysi>

(where ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX is defined in the file
Trilinos/cmake/std/atdm/<system_name>/all_supported_builds.sh) and the full
driver script name is:

  Trilinos/cmake/ctest/drivers/atdm/<system_name>/drivers/
    ${ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX}<build-name-keysi>.sh

For example, ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX=Trilinos-atdm-ats2- and
<build-name-keysi>=cuda-9.2-release-debug gives the full build name
'Trilinos-atdm-ats2-cuda-9.2-release-debug'.

However, if env var:

  ATDM_CTEST_S_USE_FULL_BUILD_NAME=1

is set, then <build-name-keysi> must match the full build name
(e.g. <build-name-keysi>=Trilinos-atdm-ats2-9.2_fpic_opt) giving the
full script name:

  Trilinos/cmake/ctest/drivers/atdm/<system_name>/drivers/
    <build-name-keysi>.sh

Tail -f the generated files <full_build_name>/smart-jenkins-driver.out to see
details of each run.

To select the default env to load instead of 'default', use:

  env ATDM_CTEST_S_DEFAULT_ENV=<system_name>-default \\
  ./ctest-s-local-test-driver.sh <build-name-1> >build-name-2> ...

(For example, one must set ATDM_CTEST_S_DEFAULT_ENV=cee-rhel7-default to run
the 'cee-rhel7' builds on CEE RHE7 machines. Otherwise the 'sems-rhel7' env
will be selected which is the default env on those machines.)

To control the list of packages tested, to build from scratch (default is to
rebuild), and not submit to CDash, use, for example:

  env \\
    Trilinos_PACKAGES=<pkg0>,<pkg1>,... \\
    CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=TRUE \\
    CTEST_DO_SUBMIT=OFF \\
  ./ctest-s-local-test-driver.sh <build-name-1> <build-name-2> ...

To test local installs, one can also set env vars:

  ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX=install
  CTEST_DO_INSTALL=ON

That will cause Trilinos to be installed into a new install/ directory under
the build directory:

  <full-build-name>/SRC_AND_BUILD/BUILD/install/

Other options that are good to set sometimes include:

  CTEST_DO_CONFIGURE=OFF
  CTEST_DO_BUILD=OFF
  CTEST_BUILD_FLAGS='<flag0> <flag1> ...'
  CTEST_DO_TEST=OFF
  CTEST_PARALLEL_LEVEL=<N>
  CTEST_DO_SUBMIT=OFF

See the documentation for TRIBITS_CTEST_DRIVER() for more details.

Options:

  -h | --help  Print this help string
"

if [[ "$@" == "-h" ]] ||  [[ "$@" == "--help" ]]; then
  echo "$CTEST_S_LOCAL_DRIVER_HELP_STR"
  exit 0
fi


#
# Functions
#

function atdm_ctest_s_get_build_name {
  build_name_body=$1
  if [[ "${ATDM_CTEST_S_USE_FULL_BUILD_NAME}" == "1" ]] ; then
    build_name="${build_name_body}"
  else
    build_name="${ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX}${build_name_body}"
  fi
  echo ${build_name}
}


function atdm_ctest_s_assert_driver_script_exists {
  build_name_body=$1
  build_name=$(atdm_ctest_s_get_build_name ${build_name_body})
  driver_script="${ATDM_TRILINOS_DIR}/cmake/ctest/drivers/atdm/${ATDM_CONFIG_SYSTEM_NAME}/drivers/${build_name}.sh"
  if [[ ! -e "${driver_script}" ]]; then
    echo
    echo "***"
    echo "*** ERROR: The driver script:"
    echo "***"
    echo "***   ${driver_script}"
    echo "***"
    echo "*** for the specified build:"
    echo "***"
    echo "***   ${build_name_body}"
    echo "***"
    echo "*** does not exist!"
    echo "***"
    return 1
  fi
  return 0
}


function atdm_ctest_s_assert_driver_scripts_exist {
  all_driver_scripts_exist=1
  for build_name in ${ATDM_ARRAY_OF_BUILDS[@]} ; do
    atdm_ctest_s_assert_driver_script_exists ${build_name} \
      || all_driver_scripts_exist=0
  done
  if [[ "${all_driver_scripts_exist}" != "1" ]]; then
    echo
    echo "Aborting the script and not running any builds!"
    exit 1
  fi
}


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

if [ "$ATDM_CTEST_S_DEFAULT_ENV" == "" ] ; then
  ATDM_CTEST_S_DEFAULT_ENV=default
fi
#echo "ATDM_CTEST_S_DEFAULT_ENV = ${ATDM_CTEST_S_DEFAULT_ENV}"

echo
echo "Load some env to get python, cmake, etc ..."
echo
source $STD_ATDM_DIR/load-env.sh ${ATDM_CTEST_S_DEFAULT_ENV}
# NOTE: Above, it does not matter which env you load.  Any of them will
# provide the right python, cmake, etc.

#
# Get the list of builds to run
#

# Must get ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX
source $STD_ATDM_DIR/$ATDM_CONFIG_SYSTEM_NAME/all_supported_builds.sh
#echo "ATDM_CONFIG_ALL_SUPPORTED_BUILDS = '${ATDM_CONFIG_ALL_SUPPORTED_BUILDS[@]}'"

if [[ "$@" == "" ]] ; then
  echo
  echo "Error, must provide 'all' or a list of supported build names which include:"
  echo
  for build_name_body in ${ATDM_CONFIG_ALL_SUPPORTED_BUILDS[@]} ; do
    build_name=$(atdm_ctest_s_get_build_name ${build_name_body})
    echo "    ${build_name}"
  done
  echo
  echo "See --help for more details!"
  exit 1
fi

ATDM_ARRAY_OF_BUILDS=$@
if [ "${ATDM_ARRAY_OF_BUILDS}" == "all" ] ; then
  ATDM_ARRAY_OF_BUILDS=${ATDM_CONFIG_ALL_SUPPORTED_BUILDS[@]}
fi

echo
echo "Running builds:"
for build_name in ${ATDM_ARRAY_OF_BUILDS[@]} ; do
  echo "    ${build_name}"
done

atdm_ctest_s_assert_driver_scripts_exist

#
# Run the builds using the ctest -S driver script
#

BASEDIR=$PWD

ln -sf ${ATDM_TRILINOS_DIR} .

for build_name_body in ${ATDM_ARRAY_OF_BUILDS[@]} ; do

  build_name=$(atdm_ctest_s_get_build_name ${build_name_body})

  echo
  date
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

  echo
  grep "failed.*out.*of" smart-jenkins-driver.out

  cd ${BASEDIR}

done

echo
date
echo
echo "Done running all of the builds!"
