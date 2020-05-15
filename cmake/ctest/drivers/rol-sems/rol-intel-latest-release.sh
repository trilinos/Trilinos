#!/bin/bash

# This is a poor-man's driver script

DRIVER_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
echo "DRIVER_SCRIPT_DIR = '$DRIVER_SCRIPT_DIR'"

TRILINOS_DIR=`readlink -f ${DRIVER_SCRIPT_DIR}/../../../..`
echo "TRILINOS_DIR='${TRILINOS_DIR}'"

source /etc/bashrc
source $TRILINOS_DIR/cmake/load_sems_dev_env.sh "sems-intel/17.0.1" "sems-openmpi/1.10.1"

export CTEST_DASHBOARD_ROOT=$PWD

ctest -V -S $DRIVER_SCRIPT_DIR/ctest-rol-mpi-release-intel-latest.sems.cmake
