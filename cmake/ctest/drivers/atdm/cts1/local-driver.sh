#!/bin/bash -l

set +x

if [ "${SLURM_CTEST_TIMEOUT}" == "" ] ; then
  SLURM_CTEST_TIMEOUT=1:20:00
  # This is just running tests, not the entire build!
fi

if [[ "${Trilinos_ENABLE_BUILD_STATS}" == "" ]] && \
   [[ ! $JOB_NAME == *"intel"* ]] \
  ; then
  export Trilinos_ENABLE_BUILD_STATS=ON
fi
echo "Trilinos_ENABLE_BUILD_STATS='${Trilinos_ENABLE_BUILD_STATS}'"
# NOTE: That above matching is a bit fragile but it avoids needing to load a
# full env and it is good enough for driving nightly builds.  (I would never
# do this with a build name coming from a user.)

set -x

source $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver-config-build.sh

set -x

if [[ "${CTEST_DO_TEST}" != "OFF" ]] ; then
atdm_run_script_on_compute_node \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver-test.sh \
  $PWD/ctest-s-driver-test.out \
  ${SLURM_CTEST_TIMEOUT}
fi # CTEST_DO_TEST
