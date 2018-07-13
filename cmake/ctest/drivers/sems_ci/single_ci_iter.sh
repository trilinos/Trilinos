#!/bin/bash

# This is a poor-man's driver script

DRIVER_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
echo "DRIVER_SCRIPT_DIR = '$DRIVER_SCRIPT_DIR'"

TRILINOS_DIR=`readlink -f ${DRIVER_SCRIPT_DIR}/../../../..`
echo "TRILINOS_DIR='${TRILINOS_DIR}'"

source /etc/bashrc
source $TRILINOS_DIR/cmake/std/GCC-4.8.4-OpenMPI-1.10.1-MpiReleaseDebugSharedPtOpenMP_env.sh
module list

export CTEST_DASHBOARD_ROOT=$PWD

# See if there are updated files:
cd Trilinos/
./cmake/tribits/python_utils/gitdist --dist-no-color fetch
NEW_COMMITS=`./cmake/tribits/python_utils/gitdist --dist-no-color log --oneline @{u} ^HEAD | grep -v "\(Git Repo\|^$\)"`
echo "NEW_COMMITS ='$NEW_COMMITS'"
cd ..

if [ "$NEW_COMMITS" != "" ] || [ "$CI_FIRST_ITERATION" == "1" ]  ; then
  #echo ctest -V -S $DRIVER_SCRIPT_DIR/ctest_linux_mpi_debug_shared_pt_ci.sems.cmake
  ctest -V -S $DRIVER_SCRIPT_DIR/ctest_linux_mpi_debug_shared_pt_ci.sems.cmake
fi
