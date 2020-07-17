#!/bin/bash -l

set +x

# Need to load env so we define some vars
source $WORKSPACE/Trilinos/cmake/std/atdm/load-env.sh $JOB_NAME

# Make adjustments for mini build of Trilinos for SPARC
if atdm_match_buildname_keyword mini ; then
  echo "This is a mini build of Trilinos for SPARC!"
  export ATDM_CONFIG_CONFIGURE_OPTIONS_FILES=cmake/std/atdm/ATDMDevEnv.cmake,cmake/std/atdm/apps/sparc/SPARC_MiniTrilinos_PACKAGES.cmake
fi

set -x

$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
