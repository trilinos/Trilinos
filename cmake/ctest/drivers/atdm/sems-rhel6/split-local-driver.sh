#!/bin/bash -l

set -x
source $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver-config-build.sh
set -x
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver-test.sh
