# This script can be used to load the appropriate environment for the
# Python 3 Pull Request testing build on a Linux machine that has access to
# the SEMS NFS mount.

# usage: $ source PullRequestPython3TestingEnv.sh

source /projects/sems/modulefiles/utils/sems-modules-init.sh

module load sems-git/2.10.1
module load sems-gcc/7.2.0

# Load the SEMS CMake Module
# - One of the SEMS modules will load CMake 3.4.x also,
#   so this will pull in the SEMS cmake 3.10.3 version
#   for Trilinos compatibility.
module load sems-cmake/3.10.3

module load sems-ninja_fortran/1.8.2

# the one we are testing.
module load sems-python/3.5.2

# this is due to some scripts in tribits calling python2
unset PYTHONHOME
