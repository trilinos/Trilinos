# This script can be used to load the appropriate environment for the
# GCC 4.8.4 Pull Request testing build on a Linux machine that has access to
# the SEMS NFS mount.

# usage: $ source PullRequestGCC4.8.4TestingEnv.sh

# After the environment is no longer needed, it can be purged using
# $ module purge
# or Trilinos/cmake/unload_sems_dev_env.sh

source /projects/sems/modulefiles/utils/sems-archive-modules-init.sh

module load sems-archive-gcc/4.8.4
module load sems-archive-openmpi/1.10.1
module load sems-archive-git/2.10.1
module load sems-archive-boost/1.63.0/base
module load sems-archive-zlib/1.2.8/base
module load sems-archive-hdf5/1.10.6/parallel
module load sems-archive-netcdf/4.7.3/parallel
module load sems-archive-parmetis/4.0.3/parallel
module load sems-archive-scotch/6.0.3/nopthread_64bit_parallel
module load sems-archive-superlu/4.3/base

# Load the SEMS CMake Module
# - One of the SEMS modules will load CMake 3.4.x also,
#   so this will pull in the SEMS cmake 3.10.3 version
#   for Trilinos compatibility.
module load sems-archive-cmake/3.17.1
module load sems-archive-ninja_fortran/1.8.2

# we will have implicitly gotten the sems python from
# the boost module above for whaever reason - reset it
# to one that has the mock packages installed
module unload sems-archive-python
# module load sierra-python/2.7.15 - permissions do not allow this, but the execs are ok
export PATH=/projects/sierra/linux_rh7/install/Python/2.7.15/bin:${PATH}
PATH=/projects/sierra/linux_rh7/install/Python/extras/bin:${PATH}
export PYTHONPATH=/projects/sierra/linux_rh7/install/Python/extras/lib/python2.7/site-packages:${PYTHONPATH}
export MANPATH=/projects/sierra/linux_rh7/install/Python/2.7.15/share/man:${MANPATH}
unset PYTHONHOME
