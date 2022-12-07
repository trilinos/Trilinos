#!/bin/bash

. /etc/bash.bashrc.local

BUILD_JOB_ID=$(sbatch --parsable configureBuild.sh)
#run tests after configure & build has finished with return code 0
PERF_TEST_JOB_ID=$(sbatch --parsable --dependency=aftercorr:$BUILD_JOB_ID runPerfTests.sh)
sbatch --dependency=afterany:$PERF_TEST_JOB_ID pushResults.sh
