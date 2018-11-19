#!/bin/bash -l

if [ "${SRUN_CTEST_TIME_LIMIT_MINUTES}" == "" ] ; then
  SRUN_CTEST_TIME_LIMIT_MINUTES=180  # Default limit is 3 hours
fi

if [ "${Trilinos_CTEST_DO_ALL_AT_ONCE}" == "" ] ; then
  export Trilinos_CTEST_DO_ALL_AT_ONCE=TRUE
fi

set -x

srun -N 1 --constraint=k80 -J $JOB_NAME \
  --time=${SRUN_CTEST_TIME_LIMIT_MINUTES} \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh

# NOTE: 'srun' is in different places on 'hansen' and 'shiller' so you can't
# provide full path to 'srun'.  You just need to hope it is in the path.
