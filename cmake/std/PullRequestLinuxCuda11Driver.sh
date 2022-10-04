#!/bin/bash -l

if [ "${BSUB_CTEST_TIME_LIMIT}" == "" ] ; then
  export BSUB_CTEST_TIME_LIMIT=12:00
fi

if [ "${Trilinos_CTEST_DO_ALL_AT_ONCE}" == "" ] ; then
  export Trilinos_CTEST_DO_ALL_AT_ONCE=TRUE
fi

# comment out sh and add what we need individually.
#source $WORKSPACE/Trilinos/cmake/std/atdm/load-env.sh $JOB_NAME

set -x

#TODO: review appropriate job size
#cuda 10.1.105 will be built with rhel7G

bsub -x -q dev -Is -n 20 -J ${JOB_NAME} -W ${BSUB_CTEST_TIME_LIMIT} ${WORKSPACE}/Trilinos/cmake/std/PullRequestLinuxDriver.sh

#if [[ "${JOB_NAME}" == "trilinos-folder/Trilinos_pullrequest_cuda_10.1.105"* ]] ; then
#  bsub -x -Is -q rhel7G -n 16 -J ${JOB_NAME} -W ${BSUB_CTEST_TIME_LIMIT} ${WORKSPACE}/Trilinos/cmake/std/PullRequestLinuxDriver.sh
#else
#  bsub -x -Is -q rhel7F -n 16 -J ${JOB_NAME} -W ${BSUB_CTEST_TIME_LIMIT} ${WORKSPACE}/Trilinos/cmake/std/PullRequestLinuxDriver.sh
#fi
# NOTE: Above, this bsub command should grab a single rhel7F (Firestone,
# Dual-Socket POWER8, 8 cores per socket, K80 GPUs) node.  The option '-x'
# makes sure that only this job runs on that node.  The options '-n 16' and
# '-q rhel7G' should make bsub allocate a single one of these nodes.
