################################################################################
#
# Get info about a custom system confiugraiton if it is selected
#
# As input, this script looks for the following vars:
#
#   ATDM_CONFIG_BUILD_NAME
#   ATDM_CONFIG_CUSTOM_CONFIG_DIR_ARG
#   ATDM_CONFIG_REGISTER_CUSTOM_CONFIG_DIR
#
# If this script selects a custom configuration then it will set the vars:
#
#   ATDM_CONFIG_REAL_HOSTNAME
#   ATDM_CONFIG_CDASH_HOSTNAME
#   ATDM_CONFIG_SYSTEM_NAME
#   ATDM_CONFIG_SYSTEM_DIR
#   ATDM_CONFIG_CUSTOM_CONFIG_DIR
#  
# If it can't, these vars will not get set (and in fact are unset)
#
# If an error does not occur, then the var:
#
#   ATDM_CONFIG_GET_CUSTOM_SYSTEM_INFO_COMPLETED=1
#
# will get set.
#
################################################################################

# Clean out vars in case this crashes before finishing

unset ATDM_CONFIG_REAL_HOSTNAME
unset ATDM_CONFIG_CDASH_HOSTNAME
unset ATDM_CONFIG_SYSTEM_NAME
unset ATDM_CONFIG_SYSTEM_DIR
unset ATDM_CONFIG_CUSTOM_CONFIG_DIR

unset ATDM_CONFIG_GET_CUSTOM_SYSTEM_INFO_COMPLETED

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

# See if a custom confiugration directory is passed in or registered as env var
unset ATDM_CUSTOM_CONFIG_DIR
if [[ "${ATDM_CONFIG_CUSTOM_CONFIG_DIR_ARG}" != "" ]] ; then
  ATDM_CUSTOM_CONFIG_DIR=${ATDM_CONFIG_CUSTOM_CONFIG_DIR_ARG}
  ATDM_CUSTOM_CONFIG_DIR_SOURCE="second argument"
elif [[ "${ATDM_CONFIG_REGISTER_CUSTOM_CONFIG_DIR}" != "" ]]; then
  ATDM_CUSTOM_CONFIG_DIR=${ATDM_CONFIG_REGISTER_CUSTOM_CONFIG_DIR}
  ATDM_CUSTOM_CONFIG_DIR_SOURCE="ATDM_CONFIG_REGISTER_CUSTOM_CONFIG_DIR"
fi
#echo "ATDM_CUSTOM_CONFIG_DIR = '${ATDM_CUSTOM_CONFIG_DIR}'"

# Assert that the custom configuration directory exists and contains
# environment.sh and get the custom system name
unset custom_system_name
if [[ "${ATDM_CUSTOM_CONFIG_DIR}" != "" ]]; then
  if [ ! -d "${ATDM_CUSTOM_CONFIG_DIR}" ] ; then
    echo "Error, '${ATDM_CUSTOM_CONFIG_DIR}' from ${ATDM_CUSTOM_CONFIG_DIR_SOURCE} must point to a valid directory with a user-defiend configuration!"
    return
  elif [ ! -e "${ATDM_CUSTOM_CONFIG_DIR}/environment.sh" ] ; then
    echo "Error, directory '${ATDM_CUSTOM_CONFIG_DIR}' from ${ATDM_CUSTOM_CONFIG_DIR_SOURCE} exists but the file '${ATDM_CUSTOM_CONFIG_DIR}/environment.sh' does not exist!"
    return
  fi
  custom_system_name=$(basename ${ATDM_CUSTOM_CONFIG_DIR})
fi
#echo "custom_system_name = '${custom_system_name}'"

# Determine if we should select the custom configuration or not
if [[ "${custom_system_name}" == "" ]] ; then
  # No custom configuration even defined
  select_custom_system_config=0
elif [[ "${ATDM_CONFIG_CUSTOM_CONFIG_DIR_ARG}" != "" ]] ; then
  # Alway use the custom configuraiton if passed into load-env.sh
  select_custom_system_config=1
elif [[ $ATDM_CONFIG_BUILD_NAME == *"${custom_system_name}"* ]] ; then
  # Use the registered custom configuration because it was listed in the build
  # name.
  select_custom_system_config=1
else
 # Otherwise, the custom system name did not appear in the build name so don't
 # use the custom configuuration.
  select_custom_system_config=0
fi
#echo "select_custom_system_config = '${select_custom_system_config}'"

# Set varaibles for custom configuration if selected
if [[ "${select_custom_system_config}" == "1" ]] ; then
  source ${ATDM_CONFIG_SCRIPT_DIR}/utils/set_custom_system_info_vars.sh \
    ${ATDM_CUSTOM_CONFIG_DIR}
fi

export ATDM_CONFIG_GET_CUSTOM_SYSTEM_INFO_COMPLETED=1
