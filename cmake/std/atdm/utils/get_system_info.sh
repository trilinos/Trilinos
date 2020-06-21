################################################################################
#
# Get info about the selected system configuration
#
# As input, this script looks for the following env vars:
#
#   ATDM_CONFIG_BUILD_NAME
#   ATDM_CONFIG_CUSTOM_CONFIG_DIR_ARG
#   ATDM_CONFIG_REGISTER_CUSTOM_CONFIG_DIR
#
# And this script will set the env vars:
#
#   ATDM_CONFIG_REAL_HOSTNAME
#   ATDM_CONFIG_CDASH_HOSTNAME
#   ATDM_CONFIG_SYSTEM_NAME
#   ATDM_CONFIG_SYSTEM_DIR
#   ATDM_CONFIG_CUSTOM_SYSTEM_DIR  (if custom config env selected)
#
# or it will error out if it can't select a system configuration
#

# Assert this script is sourced, not run!
called=$_
if [ "$called" == "$0" ] ; then
  echo "This script '$0' is being called.  Instead, it must be sourced!"
  exit 1
fi

if [ -z "$ATDM_CONFIG_SCRIPT_DIR" ] ; then
  echo "Error, must set ATDM_CONFIG_SCRIPT_DIR in env!"
  return
fi

source ${ATDM_CONFIG_SCRIPT_DIR}/utils/unset_atdm_config_vars_system_info.sh

# First, look for a custom system configuration
source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_custom_system_info.sh
if [[ "${ATDM_CONFIG_GET_CUSTOM_SYSTEM_INFO_COMPLETED}" != "1" ]] ; then
  return
fi

# Second, try to match known system configuration
if [[ $ATDM_CONFIG_SYSTEM_NAME == "" ]] ; then
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/get_known_system_info.sh
fi

export ATDM_CONFIG_KNOWN_HOSTNAME=${ATDM_CONFIG_CDASH_HOSTNAME}  # Deprecated!
export ATDM_CONFIG_KNOWN_SYSTEM_NAME=${ATDM_CONFIG_SYSTEM_NAME}  # Deprecated!
