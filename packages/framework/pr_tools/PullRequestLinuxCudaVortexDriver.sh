#!/bin/bash -l

if [ "${BSUB_CTEST_TIME_LIMIT}" == "" ] ; then
  export BSUB_CTEST_TIME_LIMIT=12:00
fi

if [ "${Trilinos_CTEST_DO_ALL_AT_ONCE}" == "" ] ; then
  export Trilinos_CTEST_DO_ALL_AT_ONCE=TRUE
fi

# comment out sh and add what we need individually.
#source $WORKSPACE/Trilinos/packages/framework/pr_tools/atdm/load-env.sh $JOB_NAME

set -x

rdc_regex=".*(-rdc)"
if [[ ${JOB_BASE_NAME:?} =~ ${rdc_regex} ]]; then
    export TRILINOS_MAX_CORES=10
fi
bsub -Is -nnodes 2 -J ${JOB_NAME} -W ${BSUB_CTEST_TIME_LIMIT} \
  ${WORKSPACE}/Trilinos/packages/framework/pr_tools/PullRequestLinuxDriver.sh

# NOTE: Above, this bsub command should grab a two
# nodes.  The option '-x' makes sure that only this
# job runs on those nodes.
