#!/bin/bash -l

set +x

if [[ "${SALLOC_CTEST_LIMIT_MINUTES}" == "" ]] ; then
  SALLOC_CTEST_LIMIT_MINUTES=4:00:00
  # From prior builds on 'stria', it looks like 4 hours should be plenty of
  # time to do the build and tests as it has taken only about 2 1/2 hours to
  # do everything.
fi

source $WORKSPACE/Trilinos/cmake/std/atdm/load-env.sh $JOB_NAME
echo

if [[ "${ATDM_CONFIG_WCID_ACCOUNT}" == "" ]] ; then
  export ATDM_CONFIG_WCID_ACCOUNT=${ATDM_CONFIG_WCID_ACCOUNT_DEFAULT}
fi

set -x

salloc -N 1 --time=${SALLOC_CTEST_LIMIT_MINUTES} -p short,batch \
  --account=${ATDM_CONFIG_WCID_ACCOUNT} \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/van1-tx2/local-driver-on-allocation.sh

set -x

# NOTE: Above, we get a single compute-node allocation using 'salloc' and then
# run the build and tests from inside of that.

# NOTE: We might need to switch from salloc to the more complex sbatch
# appraoch in the function atdm_run_script_on_compute_node.  If we see random
# ODTE errors then that is what we will need to do.
