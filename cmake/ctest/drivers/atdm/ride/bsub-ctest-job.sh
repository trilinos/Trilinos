#! /bin/bash
#PBS -l nodes=1:ppn=64
#PBS -l cput=$BSUB_CTEST_TIME_LIMIT
#PBS -N $JOB_NAME

$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
