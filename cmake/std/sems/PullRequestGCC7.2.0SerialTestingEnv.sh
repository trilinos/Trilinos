# This script can be used to load the appropriate environment for the
# GCC 7.2.0 Pull Request testing build on a Linux machine that has access to
# the SEMS NFS mount.

# usage: $ source PullRequestGCC7.2.0TestingEnv.sh

module purge

source /projects/sems/modulefiles/utils/sems-archive-modules-init.sh

module load sems-archive-gcc/7.2.0
module load sems-archive-boost/1.63.0/base
module load sems-archive-zlib/1.2.8/base
module load sems-archive-hdf5/1.10.6/base
module load sems-archive-netcdf/4.7.3/base
module load sems-archive-metis/5.1.0/base
module load sems-archive-superlu/4.3/base

module load sems-archive-cmake/3.17.1
module load sems-archive-ninja_fortran/1.8.2

# Boost loads sems-archive-python/2.7.x as a side effect.
module unload sems-archive-python
module load sems-archive-python/3.5.2

module load sems-archive-git/2.10.1

# add the OpenMP environment variable we need
export OMP_NUM_THREADS=2
