#!/bin/bash -l

set +x

if [ "${BSUB_CTEST_TIME_LIMIT}" == "" ] ; then
  export BSUB_CTEST_TIME_LIMIT=12:00
fi

if [ "${EXCLUDE_NODES_FROM_BSUB}" == "" ] ; then
  if [ "${ATDM_CONFIG_CDASH_HOSTNAME}" == "white" ] ; then
    EXCLUDE_NODES_FROM_BSUB="-R hname!=white26&&hname!=white27"
  fi
fi

if [[ "${Trilinos_ENABLE_BUILD_STATS}" == "" ]] ; then
  export Trilinos_ENABLE_BUILD_STATS=ON
fi
echo "Trilinos_ENABLE_BUILD_STATS='${Trilinos_ENABLE_BUILD_STATS}'"

source $WORKSPACE/Trilinos/cmake/std/atdm/load-env.sh $JOB_NAME

set -x

bsub -x -Is -q $ATDM_CONFIG_QUEUE -n 16 -J $JOB_NAME -W $BSUB_CTEST_TIME_LIMIT \
  ${EXCLUDE_NODES_FROM_BSUB} \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh

# NOTE: Above, this bsub command should grab a single rhel7F (Firestone,
# Dual-Socket POWER8, 8 cores per socket, K80 GPUs) node.  The option '-x'
# makes sure that only this job runs on that node.  The options '-n 16' and
# '-q rhel7G' should make bsub allocate a single one of these nodes.
