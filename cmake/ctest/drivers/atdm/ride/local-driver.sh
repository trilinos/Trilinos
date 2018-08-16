#!/bin/bash -l

if [ "${BSUB_CTEST_TIME_LIMIT}" == "" ] ; then
  export BSUB_CTEST_TIME_LIMIT=12:00
fi

if [ "${Trilinos_CTEST_DO_ALL_AT_ONCE}" == "" ] ; then
  export Trilinos_CTEST_DO_ALL_AT_ONCE=TRUE
fi

source $WORKSPACE/Trilinos/cmake/std/atdm/load-env.sh $JOB_NAME

set -x

bsub -x -Is -q $ATDM_CONFIG_QUEUE -n 16 -J $JOB_NAME -W $BSUB_CTEST_TIME_LIMIT \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh

# NOTE: Above, this bsub command should grab a single rhel7F (Firestone,
# Dual-Socket POWER8, 8 cores per socket, K80 GPUs) node.  The option '-x'
# makes sure that only this job runs on that node.  The options '-n 16' and
# '-q rhel7G' should make bsub allocate a single one of these nodes.
