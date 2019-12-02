#
# Custom builds for van1-tx2 env
#

if   [[ $ATDM_CONFIG_BUILD_NAME == *"arm-19.2-openmpi-3.1.4"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"arm-19.2_openmpi-3.1.4"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"arm-19.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"arm-19"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"arm"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"default" ]] \
  ; then
  export ATDM_CONFIG_COMPILER=ARM-19.2_OPENMPI-3.1.4

elif [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.2.0-openmpi-3.1.4"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.2.0_openmpi-3.1.4"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.2.0"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu"* ]] \
  ; then
  export ATDM_CONFIG_COMPILER=GNU-7.2.0_OPENMPI-3.1.4

else
  echo
  echo "***"
  echo "*** ERROR: A supported compiler was not selected for 'van1-tx2' (astra) env"
  echo "***"
  echo "*** Supported compilers include:"
  echo "***"
  echo "****  arm-19.2-openmpi-3.1.4    (arm-19.2, arm, default)"
  echo "****  gnu-7.2.0-openmpi-3.1.4   (gnu-7.2.0, gnu)"
  echo "***"  
  return

fi
