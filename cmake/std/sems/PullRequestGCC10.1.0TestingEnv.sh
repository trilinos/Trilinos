# This script can be used to load the appropriate environment for the
# GCC 10.1.0 Pull Request testing build on a Linux machine that has access to
# the SEMS NFS mount.

# usage: $ source PullRequestGCC8.3.0TestingEnv.sh

module purge

source /projects/sems/modulefiles/utils/sems-modules-init.sh

module load sems-gcc/8.3.0
module load sems-git/2.11.1
module load sems-openmpi/4.0.5
module load sems-boost/1.70.0
module load sems-zlib/1.2.11
module load sems-hdf5/1.10.7
module load sems-netcdf-c/4.7.3
module load sems-netcdf-cxx/4.2
module load sems-netcdf-fortran/4.5.3
module load sems-parallel-netcdf/1.12.1
module load sems-metis-int64/5.1.0
module load sems-parmetis-int64/4.0.3
module load sems-scotch-int64/6.0.3
module load sems-superlu/4.3
module load sems-ninja/1.10.1
module load sems-cmake/3.18.4

# add the OpenMP environment variable we need
export OMP_NUM_THREADS=2
