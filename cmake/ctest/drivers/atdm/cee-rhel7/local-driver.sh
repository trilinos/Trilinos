#!/bin/bash -l

set +x

# Need to load env so we define some vars
source $WORKSPACE/Trilinos/cmake/std/atdm/load-env.sh $JOB_NAME

# Ensure that we don't set both Trilinos_PACKAGES and Trilinos_PACKAGE_ENABLES_FILE
if [ -z $Trilinos_PACKAGES ]; then
  # Make adjustments for mini build of Trilinos for SPARC
  if atdm_match_buildname_keyword mini ; then
    echo "This is a mini build of Trilinos for SPARC!"
    export ATDM_CONFIG_CONFIGURE_OPTIONS_FILES=cmake/std/atdm/apps/sparc/SPARCMiniTrilinosPackagesEnables.cmake,cmake/std/atdm/ATDMDevEnv.cmake
    export Trilinos_PACKAGE_ENABLES_FILE=$WORKSPACE/Trilinos/cmake/std/atdm/apps/sparc/SPARCMiniTrilinosPackagesEnables.cmake
    # NOTE: Above, we list SPARCMiniTrilinosPackagesEnables.cmake before
    # ATDMDevEnv.cmake so that defaults for cache vars are set there before they
    # get set in ATDMDevEnv.cmake.
  fi
fi

set -x

$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
