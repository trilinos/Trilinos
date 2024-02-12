#!/bin/bash -l

set +x

# Need to load env so we define some vars
source $WORKSPACE/Trilinos/cmake/std/atdm/load-env.sh $JOB_NAME

# Make adjustments for the XL builds
if atdm_match_buildname_keyword xl ; then
  echo "This is an XL build!"
  # For XL, do not build tests and examples by default.
  if [[ "${Trilinos_INNER_ENABLE_TESTS}" == "" ]]; then
    export Trilinos_INNER_ENABLE_TESTS=OFF
  fi
  # Don't do the cuda-aware build for the XL builds if you are not building
  # internal tests and examples.
  if [[ "${Trilinos_INNER_ENABLE_TESTS}" == "OFF" ]]; then
    export Trilinos_CTEST_RUN_GPU_AWARE_MPI=0
  fi
  # Only enable the SPARC packages by default
  export ATDM_CONFIG_CONFIGURE_OPTIONS_FILES=cmake/std/atdm/ATDMDevEnv.cmake,cmake/std/atdm/apps/sparc/SPARCTrilinosPackagesEnables.cmake

  # Ensure that we don't set both Trilinos_PACKAGES and Trilinos_PACKAGE_ENABLES_FILE
  if [ -z $Trilinos_PACKAGES ]; then
      export Trilinos_PACKAGE_ENABLES_FILE=$WORKSPACE/Trilinos/cmake/std/atdm/apps/sparc/SPARCMiniTrilinosPackagesEnables.cmake
  fi
fi

# Allow default setting for TPETRA_ASSUME_GPU_AWARE_MPI=0 in trilinos_jsrun
unset TPETRA_ASSUME_GPU_AWARE_MPI
atdm_config_ctest_regex_old="$ATDM_CONFIG_CTEST_REGEX"
export ATDM_CONFIG_CTEST_REGEX="$ATDM_CONFIG_CTEST_REGEX -E Adelus*"

echo
echo "======================================================================="
echo "***"
echo "*** Run update, configure, build and testing with non-CUDA-aware MPI"
echo "***"
echo

set -x
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ats2/local-driver-single-build.sh
set +x

export ATDM_CONFIG_CTEST_REGEX="$atdm_config_ctest_regex_old"

if [[ "${Trilinos_CTEST_RUN_GPU_AWARE_MPI}" == "1" ]]; then
  echo
  echo "======================================================================="
  echo "***"
  echo "*** Running the follow-up testing with CUDA-aware MPI"
  echo "***"
  echo
  export CTEST_BUILD_NAME=${ATDM_CONFIG_BUILD_NAME}_cuda-aware-mpi
  if [[ "${CTEST_TEST_TYPE}" == "Experimental" ]] ; then
    export CTEST_BUILD_NAME="${CTEST_BUILD_NAME}-exp"
  fi
  export TPETRA_ASSUME_GPU_AWARE_MPI=1
  export CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE
  export CTEST_DO_UPDATES=OFF
  set -x
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ats2/local-driver-single-build.sh
  set +x
fi

# NOTE: We allow the configure and build to be performed again in the
# follow-up CUDA-aware MPI build though they don't need to be.  This so that
# developers that are looking on CDash will be able to see configure-related
# information associated with that build.  We also need to run configure and
# build again to make the cdash_analyze_and_report.py tool not report this
# build as missing configure and build data.  However, we don't do an update
# again because that would always give zero updates and it would confuse
# developers looking on CDash to see what changed from yesterday's build and
# see no updates were pulled.  It is better to not do the update and not
# upload the empty Update.txt file and avoid that confusion.

# NOTE: Above, we put an '-exp' on the end of the cuda-aware-mpi build name to
# make it more clear that this is an experimental build!
