#!/bin/bash

export WATCHR_BUILD_NAME="Vortex CUDA"
source $TRILINOS_SRC/cmake/std/atdm/load-env.sh cuda-10.1.243_gnu-7.3.1_spmpi-2019.06.24-release
cd $TRILINOS_SRC
export TRILINOS_GIT_SHA=`git rev-parse HEAD`
export TPETRA_ASSUME_CUDA_AWARE_MPI=1

cd $WORKSPACE/build

#Don't fail the whole Jenkins build if tests fail. There will just
#be a gap in the data series for failing tests.
ctest -V || true
