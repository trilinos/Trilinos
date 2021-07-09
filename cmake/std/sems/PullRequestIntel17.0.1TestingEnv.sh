# This script can be used to load the appropriate environment for the
# GCC 4.9.3 Pull Request testing build on a Linux machine that has access to
# the SEMS NFS mount.

# usage: $ source PullRequestIntel17.0.1TestingEnv.sh

# After the environment is no longer needed, it can be purged using
# $ module purge
# or Trilinos/cmake/unload_sems_dev_env.sh

module purge

source /projects/sems/modulefiles/utils/sems-modules-init.sh

export SEMS_FORCE_LOCAL_COMPILER_VERSION=4.9.3
module load sems-gcc/5.3.0
module load sems-intel/17.0.1
module load sems-mpich/3.2
module load sems-boost/1.63.0/base
module load sems-zlib/1.2.8/base
module load sems-hdf5/1.10.6/parallel
module load sems-netcdf/4.7.3/parallel
module load sems-parmetis/4.0.3/parallel
module load sems-scotch/6.0.3/nopthread_64bit_parallel
module load sems-superlu_dist/5.2.2/32bit_parallel

module load sems-cmake/3.17.1
module load sems-ninja_fortran/1.8.2

module load sems-git/2.10.1

module unload sems-python
module load sems-python/3.5.2

# add the OpenMP environment variable we need
export OMP_NUM_THREADS=2

# required for cmake > 3.10 during configure compiler testing.
export LDFLAGS=-lifcore
