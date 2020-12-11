################################################################################
#
# Source to set up the env to do ATDM configuration of Trilinos.
#
################################################################################

# Assert this script is sourced, not run!
if [ "${BASH_SOURCE[0]}" == "${0}" ] ; then
  echo "This script '${0}' is being called.  Instead, it must be sourced!"
  exit 1
fi

# Return the absoute directory of some relative directory path.
#
# This uses a temp shell to cd into the directory and then uses pwd to get the
# path.
function atdm_config_get_abs_dir_path() {
  [ -z "$1" ] && { pwd; return; }
  (cd -P -- "$1" && pwd)
}
export atdm_config_get_abs_dir_path

# Get the base dir for the sourced script
ATDM_SCRIPT_DIR=$(readlink -f \
                $(echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"))
#echo "ATDM_SCRIPT_DIR = '$ATDM_SCRIPT_DIR'"

# Absolute path to this scripts dir
export ATDM_CONFIG_SCRIPT_DIR=$ATDM_SCRIPT_DIR


#
# A) Read the command-line arguments
#

# Make sure job-name is passed in as first (and only ) arguemnt
if [ "$1" == "" ] ; then
  echo "Error, the first argument must be the build name with keyword names!"
  return
fi

# Get the build name as 1st argument
export ATDM_CONFIG_BUILD_NAME=$1
export ATDM_CONFIG_JOB_NAME=$ATDM_CONFIG_BUILD_NAME  # Deprecated!

# Optional user-defined system configuration dir as 2nd argument
unset ATDM_CONFIG_CUSTOM_CONFIG_DIR_ARG
if [ "$2" != "" ] ; then
  export ATDM_CONFIG_CUSTOM_CONFIG_DIR_ARG=$2
fi

#
# B) Get the host name, system name, and system configuration diecrory
#

source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_system_info.sh

if [[ $ATDM_CONFIG_SYSTEM_NAME == "" ]] ; then
  echo "Error, could not determine a system configuration for hostname='$realHostname', aborting env loading script!"
  return
fi

#
# C) Set Trilinos base directory
#

# Get the Trilins base dir
export ATDM_CONFIG_TRILINOS_DIR=`atdm_config_get_abs_dir_path ${ATDM_CONFIG_SCRIPT_DIR}/../../..`
if [[ $ATDM_CONFIG_VERBOSE == "1" ]] ; then
  echo "ATDM_CONFIG_TRILINOS_DIR = $ATDM_CONFIG_TRILINOS_DIR"
fi

#
# D) Parse $ATDM_CONFIG_BUILD_NAME for consumption by the system-specific
# environoment.sh script
#

source ${ATDM_CONFIG_SCRIPT_DIR}/utils/unset_atdm_config_vars_build_options.sh

source ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_build_options.sh

#
# E) Load the matching env
#

# Set other vaues to empty by default
source ${ATDM_CONFIG_SCRIPT_DIR}/utils/unset_atdm_config_vars_environment.sh

# Abort if the set_build_options script aborted (but after enviornment.sh vars unset)!
if [ "${ATDM_CONFIG_FINISHED_SET_BUILD_OPTIONS}" != "1" ]; then
  return
fi

# Set the location for NVCC wrapper for source dir unless this is from an
# install of Trilinos!
if [[ "${ATDM_SCRIPT_DIR}" == *"atdm-trilinos" ]] ; then
  # This is being invoked from an install so use the installed nvcc_wrapper we
  # installed in this installation directory
  export ATDM_CONFIG_NVCC_WRAPPER="${ATDM_SCRIPT_DIR}/nvcc_wrapper"
else
  # This is from the Trilinos source tree so grab nvcc_wrapper from there!
  export ATDM_CONFIG_NVCC_WRAPPER="${ATDM_CONFIG_TRILINOS_DIR}/packages/kokkos/bin/nvcc_wrapper"
fi

source ${ATDM_CONFIG_SCRIPT_DIR}/utils/atdm_config_helper_funcs.sh

source ${ATDM_CONFIG_SYSTEM_DIR}/environment.sh

if [ "$ATDM_CONFIG_COMPLETED_ENV_SETUP" != "TRUE" ] ; then
  echo
  echo "***"
  echo "*** ERROR: Environment setup was not successful, see above errors!"
  echo "***"
  echo
  return
fi

#
# F) Override parallel build and ctest levels
#

if [ "${ATDM_CONFIG_BUILD_COUNT_OVERRIDE}" != "" ] ; then
  export ATDM_CONFIG_BUILD_COUNT=${ATDM_CONFIG_BUILD_COUNT_OVERRIDE}
fi

if [ "${ATDM_CONFIG_PARALLEL_COMPILE_JOBS_LIMIT_OVERRIDE}" != "" ] ; then
  export ATDM_CONFIG_PARALLEL_COMPILE_JOBS_LIMIT=${ATDM_CONFIG_PARALLEL_COMPILE_JOBS_LIMIT_OVERRIDE}
fi

if [ "${ATDM_CONFIG_PARALLEL_LINK_JOBS_LIMIT_OVERRIDE}" != "" ] ; then
  export ATDM_CONFIG_PARALLEL_LINK_JOBS_LIMIT=${ATDM_CONFIG_PARALLEL_LINK_JOBS_LIMIT_OVERRIDE}
fi

if [ "${ATDM_CONFIG_CTEST_PARALLEL_LEVEL_OVERRIDE}" != "" ] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=${ATDM_CONFIG_CTEST_PARALLEL_LEVEL_OVERRIDE}
fi

#
# G) Set install-related stuff
#

if  [[ "${ATDM_CONFIG_USE_INSTALL_PBP_RUNNER_DEFAULT}" == "" ]] \
  && [[ "${ATDM_CONFIG_USE_JENKINS_INSTALL_DEFAULTS}" == "1" ]] ; then
  export ATDM_CONFIG_USE_INSTALL_PBP_RUNNER_DEFAULT=1
fi

if [[ "${ATDM_CONFIG_INSTALL_PBP_RUNNER}" == "" ]] \
  && [[ "${ATDM_CONFIG_INSTALL_PBP_RUNNER_DEFAULT}" != "" ]] \
  && [[ "${ATDM_CONFIG_USE_INSTALL_PBP_RUNNER_DEFAULT}" == "1" ]] ; then
  export ATDM_CONFIG_INSTALL_PBP_RUNNER="${ATDM_CONFIG_INSTALL_PBP_RUNNER_DEFAULT}"
fi

# NOTE: The ATDMDevEnvSettings.cmake module when processed will assert that
# many of these vars are correctly set.
