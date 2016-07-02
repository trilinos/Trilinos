#!/bin/bash

# Test (and push) Trilinos on any machine that has the SEMS Dev Env available
#
# To use, just symlink this script into the directory where you want to run it with:
#
#   $ cd <some-build-dir>/
#   $ ln -s <trilinos-src-dir>/cmake/std/sems/checkin-test-sems.sh .
#
# For local builds, run with:
#
#   $ ./checkin-test-sems.sh [other options] --local-do-all
#
# For testing and pushing, run with:
#
#   $ ./checkin-test-sems.sh [other options] --do-all --push
#
# NOTE: After the first time running this script (e.g. ./checkin-test-sems.sh
# --help), please modify the created file:
#
#    ./local-checkin-test-defaults.py
#
# to adjust the defaults for usage on your machine.  For example, you will
# want to adjust the parallel level -j<N> option to take better advantage of
# the cores on your local machine.  Also, youe will want to make
# --ctest-timeout larger if the machine you are on is particularly slow.
#

#
# A) Get the path to the Trilinos source directory
#
# NOTE: Script must be symlinked into the Trilinos source tree for this to
# work correctly.
#

if [ "$TRILINOS_DIR" == "" ] ; then
  _ABS_FILE_PATH=`readlink -f $0`
  _SCRIPT_DIR=`dirname $_ABS_FILE_PATH`
  TRILINOS_DIR=$_SCRIPT_DIR/../../..
fi

TRILINOS_DIR_ABS=$(readlink -f $TRILINOS_DIR)
echo "TRILINOS_DIR_ABS = $TRILINOS_DIR_ABS"

# Make sure the right env is loaded!
export TRILINOS_SEMS_DEV_ENV_VERBOSE=1
source $TRILINOS_DIR_ABS/cmake/load_ci_sems_dev_env.sh

#
# B) Set up the bulid configurations
#

echo "
-DTrilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON
-DTrilinos_TRACE_ADD_TEST=ON
" > COMMON.config
# Above, we want to allow disabling enabled downstream packages if a
# dependency turns them off.  Also, we want to see the tracing of the tests
# being added to help debug why a given test is being added or not.

# NOTE: The PT --default-builds MPI_DEBUG and SERIAL_RELEASE are already set
# up automatically!

echo "
-DTPL_ENABLE_MPI=ON
-DCMAKE_BUILD_TYPE=RELEASE
-DTrilinos_ENABLE_DEBUG=ON
-DBUILD_SHARED_LIBS=ON
-DIfpack2_Cheby_belos_MPI_1_DISABLE=ON
" > MPI_RELEASE_DEBUG_SHARED_ST.config

echo "
-DTPL_ENABLE_MPI=OFF
-DCMAKE_BUILD_TYPE=RELEASE
-DBUILD_SHARED_LIBS=ON
" > SERIAL_RELEASE_SHARED_ST.config

echo "
-DTPL_ENABLE_MPI=ON
-DCMAKE_BUILD_TYPE=RELEASE
-DTrilinos_ENABLE_DEBUG=ON
-DBUILD_SHARED_LIBS=OFF
-DTPL_FIND_SHARED_LIBS=OFF
-DTrilinos_LINK_SEARCH_START_STATIC=ON
-DIfpack2_Cheby_belos_MPI_1_DISABLE=ON
" > MPI_RELEASE_DEBUG_STATIC_ST.config

echo "
-DTPL_ENABLE_MPI=OFF
-DCMAKE_BUILD_TYPE=RELEASE
-DBUILD_SHARED_LIBS=OFF
-DTPL_FIND_SHARED_LIBS=OFF
-DTrilinos_LINK_SEARCH_START_STATIC=ON
" > SERIAL_RELEASE_STATIC_ST.config

#
# C) Create a default local defaults file
#
#  After this is created, then modify the created file for your preferences!
# 

_LOCAL_CHECKIN_TEST_DEFAULTS=local-checkin-test-defaults.py
if [ -f $_LOCAL_CHECKIN_TEST_DEFAULTS ] ; then
  echo "File $_LOCAL_CHECKIN_TEST_DEFAULTS already exists, leaving it!"
else
  echo "Creating default file $_LOCAL_CHECKIN_TEST_DEFAULTS!"
  echo "
defaults = [
  \"-j4\",
  \"--ctest-timeout=180\",
  \"--disable-packages=PyTrilinos,Pliris,Claps,TriKota\",
  \"--skip-case-no-email\",
  ]
  " > $_LOCAL_CHECKIN_TEST_DEFAULTS
fi
# Note the above defaults:
#
# -j4: Uses very few processes by default
#
# --ctest-timeout=180: A default 3 minute timeout
#
# --disable-packages: These a packages that most Trilinos developers wil not
#    want to test by default.
#
# --skip-case-no-email: Don't send email for skipped test cases
#
# To exclude tests, include the option:
#
#  \"--ctest-options=-E '(<test1>|<test2>|...)\",

#
# D) Run the checkin-test.py script!
#
# NOTE: default args are read in from the local-checkin-test-defaults.py file!
#

$TRILINOS_DIR/checkin-test.py \
"$@"
