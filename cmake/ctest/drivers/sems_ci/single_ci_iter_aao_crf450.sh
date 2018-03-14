#!/bin/bash

# This is a poor-man's driver script for all-at-once CI build on crf450 (see
# below!)

DRIVER_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
echo "DRIVER_SCRIPT_DIR = '$DRIVER_SCRIPT_DIR'"

TRILINOS_DIR=`readlink -f ${DRIVER_SCRIPT_DIR}/../../../..`
echo "TRILINOS_DIR='${TRILINOS_DIR}'"

source /etc/bashrc
TRILINOS_SEMS_DEV_ENV_VERBOSE=1
source $TRILINOS_DIR/cmake/load_sems_dev_env.sh
export PATH=/home/vera_env/common_tools/cmake-master-20170917-214d0ce/bin:$PATH

export CTEST_DASHBOARD_ROOT=$PWD

env http_proxy= \
ctest -V -S $DRIVER_SCRIPT_DIR/ctest_linux_mpi_debug_shared_pt_ci_aao.sems.cmake

# ToDo: Put back the simple CI logic based on changed file!
