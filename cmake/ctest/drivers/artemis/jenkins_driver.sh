#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on typhon: `date`"
echo

# snapshot_kokkos_into_trilinos <KOKKOS-BRANCH>
snapshot_kokkos_into_trilinos() {
    cd $WORKSPACE/Trilinos && git reset --hard HEAD && cd -
    cd $WORKSPACE/kokkos && git reset --hard origin/$1 && cd -
    $WORKSPACE/kokkos/config/snapshot.py -n $WORKSPACE/kokkos $WORKSPACE/Trilinos/packages || { echo "SNAPSHOT FAILED!" && exit 1; }
    export KOKKOS_BRANCH=$1
}

#
# TrilinosDriver settings:
#

export TDD_PARALLEL_LEVEL=2

# Trilinos settings:
#

# Submission mode for the *TrilinosDriver* dashboard
export TDD_CTEST_TEST_TYPE=Nightly
export TRIBITS_TDD_USE_SYSTEM_CTEST=1

#export CTEST_DO_SUBMIT=FALSE
#export CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE

# Machine specific environment
#

source /projects/modulefiles/utils/sems-modules-init.sh
source /projects/modulefiles/utils/kokkos-modules-init.sh

module load python/2.7.9
module load cuda/6.5.14
module load cmake/2.8.11
module load git

export FROM_JENKINS=1
export TDD_HTTP_PROXY="http://sonproxy.sandia.gov:80"
export http_proxy="http://sonproxy.sandia.gov:80"
export CUDA_LAUNCH_BLOCKING=1
export OMP_NUM_THREADS=2

# Machine independent cron_driver:
#
#openmpi-1.7-cuda6

SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`

module load intel/15.0.2/openmpi/1.8.7/cuda/6.5.14
module load superlu/4.3/intel/15.0.2/base

export KOKKOS_BRANCH=master
$SCRIPT_DIR/../cron_driver.py

snapshot_kokkos_into_trilinos develop
$SCRIPT_DIR/../cron_driver.py

# module unload intel/15.0.2/openmpi/1.8.7/cuda/6.5.14
# module unload superlu/4.3/intel/15.0.2/base

# module load gcc/4.8.4/openmpi/1.8.7/cuda/6.5.14
# module load superlu/4.3/gcc/4.8.4/base

# cd $WORKSPACE/Trilinos && git reset --hard HEAD && cd -
# $SCRIPT_DIR/../cron_driver.py

# snapshot_kokkos_into_trilinos develop
# $SCRIPT_DIR/../cron_driver.py

git status

cd $WORKSPACE/Trilinos && git reset --hard HEAD && cd -

echo
echo "Ending nightly Trilinos development testing on typhon: `date`"
echo
