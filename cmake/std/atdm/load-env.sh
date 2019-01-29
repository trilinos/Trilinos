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

# Absolute path to this scripts dir
export ATDM_CONFIG_SCRIPT_DIR=`readlink -f ${ATDM_SCRIPT_DIR}`

#
# A) Parse the command-line arguments
#

# Make sure job-name is passed in as first (and only ) arguemnt
if [ "$1" == "" ] ; then
  echo "Error, the first argument must be the build name with keyword names!"
  return
fi

# Make sure there are no other command-line arguments set
if [ "$2" != "" ] ; then
  echo "Error, this source script only accepts a single comamnd-line argument!"
  return
fi

export ATDM_CONFIG_BUILD_NAME=$1

# Set old name for backward compatiblity
export ATDM_CONFIG_JOB_NAME=$ATDM_CONFIG_BUILD_NAME

#
# B) Get the system name from the hostname
#

source ${ATDM_CONFIG_SCRIPT_DIR}/utils/unset_atdm_config_vars_system_name.sh

source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_known_system_name.sh

if [[ $ATDM_CONFIG_KNOWN_SYSTEM_NAME == "" ]] ; then
  echo "Error, could not determine known system, aborting env loading"
  return
fi

#
# C) Set ATDM_CONFIG_BUILD_NAME and Trilinos base directory
#

# Get the Trilins base dir
export ATDM_CONFIG_TRILNOS_DIR=`get_abs_dir_path ${ATDM_CONFIG_SCRIPT_DIR}/../../..`
if [[ $ATDM_CONFIG_VERBOSE == "1" ]] ; then
  echo "ATDM_CONFIG_TRILNOS_DIR = $ATDM_CONFIG_TRILNOS_DIR"
fi

#
# D) Parse $ATDM_CONFIG_BUILD_NAME for consumption by the system-specific environoment.sh
# script
#

source ${ATDM_CONFIG_SCRIPT_DIR}/utils/unset_atdm_config_vars_build_options.sh

source ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh

#
# E) Load the matching env
#

# Set other vaues to empty by default
source ${ATDM_CONFIG_SCRIPT_DIR}/utils/unset_atdm_config_vars_environment.sh

# Set the location for NVCC wrapper for source dir unless this is from an
# install of Trilinos!
if [[ "${ATDM_SCRIPT_DIR}" == *"atdm-trilinos" ]] ; then
  # This is being invoked from an install so use the installed nvcc_wrapper we
  # installed in this installation directory
  export ATDM_CONFIG_NVCC_WRAPPER="${ATDM_SCRIPT_DIR}/nvcc_wrapper"
else
  # This is from the Trilinos source tree so grab nvcc_wrapper from there!
  export ATDM_CONFIG_NVCC_WRAPPER="${ATDM_CONFIG_TRILNOS_DIR}/packages/kokkos/bin/nvcc_wrapper"
fi

source ${ATDM_CONFIG_SCRIPT_DIR}/utils/atdm_config_helper_funcs.sh

source ${ATDM_CONFIG_SCRIPT_DIR}/$ATDM_CONFIG_KNOWN_SYSTEM_NAME/environment.sh

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
