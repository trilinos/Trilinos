#!/bin/bash -l

if [ "${BSUB_CTEST_TIME_LIMIT}" == "" ] ; then
  export BSUB_CTEST_TIME_LIMIT=04:00:00
fi

# We have to generate the file bsub-ctest-job.sh because it has values in
# comments that we need to produce.  Just a regular bash script will not
# evaluate $<var_name> in a comment.

echo "#! /bin/bash
#PBS -l nodes=1:ppn=64
#PBS -l cput=$BSUB_CTEST_TIME_LIMIT
#PBS -N $JOB_NAME

$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
" > bsub-ctest-job.sh

chmod a+x bsub-ctest-job.sh

set -x

bsub -x -I $PWD/bsub-ctest-job.sh
