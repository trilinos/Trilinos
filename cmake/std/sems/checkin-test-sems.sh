#!/bin/bash

# Used to test Trilinos on any machine that has the SEMS Dev Env available.

if [ "$TRILINOS_DIR" == "" ] ; then
  _ABS_FILE_PATH=`readlink -f $0`
  _SCRIPT_DIR=`dirname $_ABS_FILE_PATH`
  TRILINOS_DIR=$_SCRIPT_DIR/../../..
fi

TRILINOS_DIR_ABS=$(readlink -f $TRILINOS_DIR)
echo "TRILINOS_DIR_ABS = $TRILINOS_DIR_ABS"

# Packages in Trilinos to disable (mostly for auotmated CI server)
DISABLE_PACKAGES=PyTrilinos,Pliris,Claps,TriKota

# Make sure the right env is loaded!
export TRILINOS_SEMS_DEV_ENV_VERBOSE=1
source $TRILINOS_DIR_ABS/cmake/load_ci_sems_dev_env.sh

echo "
-DTrilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON
-DTrilinos_TRACE_ADD_TEST=ON
" > COMMON.config

#
# Built-in Primary Tested (PT) --default-builds (DO NOT MODIFY)
#

#
# Standard Secondary Tested (ST) --st-extra-builds (DO NOT MODIFY)
#

echo "
-DCMAKE_BUILD_TYPE=RELEASE
-DTrilinos_ENABLE_DEBUG=ON
-DTPL_ENABLE_MPI=ON
" > MPI_RELEASE_DEBUG_ST.config

echo "
-DCMAKE_BUILD_TYPE=RELEASE
-DTrilinos_ENABLE_DEBUG=OFF
-DTPL_ENABLE_MPI=OFF
" > SERIAL_RELEASE_ST.config

#
# Invocation
#

$TRILINOS_DIR/checkin-test.py \
--st-extra-builds=MPI_RELEASE_DEBUG_ST,SERIAL_RELEASE_ST.config \
"$@"
