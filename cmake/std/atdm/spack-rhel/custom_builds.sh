#
# Custom builds for spack-rhel env
#
# NOTE: This file gets sourced in atdm/utils/set_build_options.sh before the
# default grep logic is applied.
#

if   [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.2.0-openmpi-1.10.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.2.0_openmpi-1.10.2"* ]] \
  ; then
  export ATDM_CONFIG_COMPILER=GNU-7.2.0_OPENMPI-1.10.2

elif [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.2.0-openmpi-2.1.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.2.0_openmpi-2.1.2"* ]] \
  ; then
  export ATDM_CONFIG_COMPILER=GNU-7.2.0_OPENMPI-2.1.2

elif [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.2.0"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"default" ]] \
  ; then
  export ATDM_CONFIG_COMPILER=GNU-7.2.0_OPENMPI-1.10.2

else
  echo
  echo "***"
  echo "*** ERROR: A supported compiler was not selected for 'spack-rhel' env"
  echo "***"
  echo "*** Supported compilers include:"
  echo "***"
  echo "****  gnu-7.2.0-openmpi-1.10.2 (default)"
  echo "****  gnu-7.2.0-openmpi-2.1.2"
  echo "***"  
  return

fi
