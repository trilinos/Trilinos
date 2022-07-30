#!/bin/bash
#SBATCH -A cfd116
#SBATCH -J trilinos-perf-tests
#SBATCH -t 01:00:00
#SBATCH -N 1
#SBATCH --mail-type=BEGIN,END,INVALID_DEPEND

TRILINOS_SRC="$HOME/spock/sources/trilinos/Trilinos-muelu"
TRILINOS_BUILD="$HOME/spock/builds/performance"

export WATCHR_BUILD_NAME="SPOCK HIP"
export WATCHR_PERF_DIR="$HOME/spock/trilinos-kernel-performance/SpockData"

source "${TRILINOS_SRC}/packages/tpetra/scripts/PerformanceTesting/spock/load_modules.sh"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

cd $TRILINOS_SRC
# watchr can embed SHA with timing data
export TRILINOS_GIT_SHA=`git rev-parse HEAD`

cd $TRILINOS_BUILD

#Don't fail the whole Jenkins build if tests fail. There will just
#be a gap in the data series for failing tests.
ctest -V
