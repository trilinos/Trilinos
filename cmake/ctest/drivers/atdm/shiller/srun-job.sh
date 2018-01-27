#!/bin/bash

if [ "${WORKSPACE}" == ""  ] ; then
  echo "Error, must set WORKSPACE var before calling"
  exit 1
fi

# Shiller/hansen settings
export SLURM_TASKS_PER_NODE=32
export CUDA_LAUNCH_BLOCKING=1
export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1
export OMP_NUM_THREADS=2
ulimit -c 0

export SUBDIR=SRC_AND_BUILD
if [ ! -e $SUBDIR ] ; then
  echo "Making $SUBDIR"
  mkdir $SUBDIR
fi

cd $SUBDIR/
echo "Current dir: $PWD"

source $WORKSPACE/Trilinos/cmake/std/atdm/shiller/environment.sh
echo
module list

ctest -V -S \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/shiller/ctest-driver.cmake
