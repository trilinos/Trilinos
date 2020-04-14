#
# Custom builds for van1-tx2 env
#

if   [[ $ATDM_CONFIG_BUILD_NAME == *"arm-20.0-openmpi-4.0.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"arm-20.0_openmpi-4.0.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"arm-20.0"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"arm-20"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"default" ]] \
  ; then
  export ATDM_CONFIG_COMPILER=ARM-20.0_OPENMPI-4.0.2

else
  echo
  echo "***"
  echo "*** ERROR: A supported compiler was not selected for 'van1-tx2' (stria) env"
  echo "***"
  echo "*** Supported compilers include:"
  echo "***"
  echo "****  arm-20.0-openmpi-4.0.2    (arm-20.0, default)"
  echo "***"
  return

fi
