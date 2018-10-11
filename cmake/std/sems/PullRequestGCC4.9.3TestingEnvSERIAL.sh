# This script can be used to load the appropriate environment for the
# GCC 4.9.3 Pull Request testing build on a Linux machine that has access to
# the SEMS NFS mount.

# usage: $ source PullRequestGCC4.9.3TestingEnv.sh

source /projects/sems/modulefiles/utils/sems-modules-init.sh

module load sems-gcc/4.9.3
module load sems-python/2.7.9
module load sems-git/2.10.1
module load sems-boost/1.63.0/base
module load sems-zlib/1.2.8/base
module load sems-hdf5/1.8.12/base
module load sems-netcdf/4.4.1/exo
module load sems-metis/5.1.0/base
# module load sems-scotch/6.0.3/nopthread_64bit_parallel
module load sems-superlu/4.3/base

# Load the SEMS CMake Module
# - One of the SEMS modules will load CMake 3.4.x also,
#   so this will pull in the SEMS cmake 3.10.3 version
#   for Trilinos compatibility.
module load sems-cmake/3.10.3

# Using CMake and Ninja modules from the ATDM project space.
# SEMS does not yet supply a recent enough version of CMake
# for the single configure/build/test capability. We are also
# using a custom version of Ninja (with Fortran support not
# available in main-line Ninja) to significantly speed up
# compile and link times.
module load atdm-env
module load atdm-cmake/3.11.1
module load atdm-ninja_fortran/1.7.2


