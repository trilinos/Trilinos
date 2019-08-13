################################################################################
#
# Get info about a custom env directory
#
################################################################################

# Clean out vars in case this crashes before finishing

unset ATDM_CONFIG_REAL_HOSTNAME
unset ATDM_CONFIG_KNOWN_HOSTNAME
unset ATDM_CONFIG_KNOWN_SYSTEM_NAME

# Assert this script is sourced, not run!
called=$_
if [ "$called" == "$0" ] ; then
  echo "This script '$0' is being called.  Instead, it must be sourced!"
  exit 1
fi

# Assert that ATDM_CONFIG_BUILD_NAME is set!
if [ -z "$ATDM_CONFIG_BUILD_NAME" ] ; then
  echo "Error, must set ATDM_CONFIG_BUILD_NAME in env!"
  return
fi

# Assert that ATDM_CONFIG_CUSTOM_CONFIG_DIR_PATH is set!
if [ -z "$ATDM_CONFIG_CUSTOM_CONFIG_DIR_PATH" ] ; then
  echo "Error, must set ATDM_CONFIG_CUSTOM_CONFIG_DIR_PATH in env!"
  return
fi

export ATDM_CONFIG_REAL_HOSTNAME=`hostname`

custom_system_name=$(basename ${ATDM_CONFIG_CUSTOM_CONFIG_DIR_PATH})

export ATDM_CONFIG_KNOWN_HOSTNAME=$ATDM_HOSTNAME
export ATDM_CONFIG_KNOWN_SYSTEM_NAME=$custom_system_name
