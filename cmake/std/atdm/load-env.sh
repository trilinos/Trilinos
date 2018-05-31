################################################################################
#
# Source to set up the env to do ATDM configuration of Trilinos.
#
################################################################################

# Assert this script is sourced, not run!
called=$_
if [ "$called" == "$0" ] ; then
  echo "This script '$0' is being called.  Instead, it must be sourced!"
  exit 1
fi

# Return the absoute directory of some relative directory path.
#
# This uses a temp shell to cd into the directory and then uses pwd to get the
# path.
function get_abs_dir_path() {
  [ -z "$1" ] && { pwd; return; }
  (cd -P -- "$1" && pwd)
}

# Get the base dir for the sourced script
ATDM_SCRIPT_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
#echo "ATDM_SCRIPT_DIR = '$ATDM_SCRIPT_DIR'"

#
# A) Parse the command-line arguments
#

# Make sure job-name is passed in as first (and only ) arguemnt
if [ "$1" == "" ] ; then
  echo "Error, the first argument must be the job name with keyword names!"
  return
fi

# Make sure there are no other command-line arguments set
if [ "$2" != "" ] ; then
  echo "Error, this source script only accepts a single comamnd-line argument!"
  return
fi

#
# B) Get the system name from the hostname
#

source $ATDM_SCRIPT_DIR/utils/get_known_system_name.sh

if [[ $ATDM_CONFIG_KNOWN_SYSTEM_NAME == "" ]] ; then
  echo "Error, could not determine known system, aborting env loading"
  return
fi

#
# C) Set JOB_NAME and Trilinos base directory
#

export JOB_NAME=$1

# Get the Trilins base dir
export ATDM_CONFIG_TRILNOS_DIR=`get_abs_dir_path $ATDM_SCRIPT_DIR/../../..`
echo "ATDM_CONFIG_TRILNOS_DIR = $ATDM_CONFIG_TRILNOS_DIR"

#
# D) Parse $JOB_NAME for consumption by the system-specific environoment.sh
# script
#

source $ATDM_SCRIPT_DIR/utils/set_build_options.sh

#
# E) Load the matching env
#

# Set other vaues to empty by default
unset OMP_NUM_THREADS
unset OMP_PROC_BIND
unset OMP_PLACES
unset OMPI_CC
unset OMPI_CXX
unset OMPI_FC
unset ATDM_CONFIG_USE_NINJA
unset ATDM_CONFIG_BUILD_COUNT
unset ATDM_CONFIG_CMAKE_JOB_POOL_LINK
unset ATDM_CONFIG_CTEST_PARALLEL_LEVEL
unset ATDM_CONFIG_KOKKOS_ARCH
unset ATDM_CONFIG_BLAS_LIB
unset ATDM_CONFIG_LAPACK_LIB
unset ATDM_CONFIG_USE_HWLOC
unset ATDM_CONFIG_HWLOC_LIBS
unset ATDM_CONFIG_USE_CUDA
unset ATDM_CONFIG_HDF5_LIBS
unset ATDM_CONFIG_NETCDF_LIBS
unset ATDM_CONFIG_MPI_EXEC
unset ATDM_CONFIG_MPI_PRE_FLAGS
unset ATDM_CONFIG_MPI_POST_FLAG
unset ATDM_CONFIG_COMPLETED_ENV_SETUP

source $ATDM_SCRIPT_DIR/$ATDM_CONFIG_KNOWN_SYSTEM_NAME/environment.sh

if [ "${ATDM_CONFIG_BUILD_COUNT_OVERRIDE}" != "" ] ; then
  export ATDM_CONFIG_BUILD_COUNT=${ATDM_CONFIG_BUILD_COUNT_OVERRIDE}
fi

if [ "${ATDM_CONFIG_CMAKE_JOB_POOL_LINK_OVERRIDE}" != "" ] ; then
  export ATDM_CONFIG_CMAKE_JOB_POOL_LINK=${ATDM_CONFIG_CMAKE_JOB_POOL_LINK_OVERRIDE}
fi

if [ "${ATDM_CONFIG_CTEST_PARALLEL_LEVEL_OVERRIDE}" != "" ] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=${ATDM_CONFIG_CTEST_PARALLEL_LEVEL_OVERRIDE}
fi

if [ "$ATDM_CONFIG_COMPLETED_ENV_SETUP" != "TRUE" ] ; then
  echo
  echo "***"
  echo "*** ERROR: Environment setup was not successful, see above errors!"
  echo "***"
  echo
  return
fi

# NOTE: The ATDMDevEnv.cmake module when processed will assert that all of
# these are set!
