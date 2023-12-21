#!/bin/bash -e

export WATCHR_BUILD_NAME="Vortex CUDA"
source $WORKSPACE/PerfScripts/loadEnv.sh
cd $TRILINOS_SRC
export TRILINOS_GIT_SHA=`git rev-parse HEAD`

cd $WORKSPACE/build

#Don't fail the whole Jenkins build if tests fail. There will just
#be a gap in the data series for failing tests.
ctest -V || true

