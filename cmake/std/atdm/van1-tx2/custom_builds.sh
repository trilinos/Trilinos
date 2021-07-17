#
# Custom builds for van1-tx2 env
#

if atdm_match_any_buildname_keyword \
     arm-20.1-openmpi-4.0.5 \
     arm-20.1_openmpi-4.0.5 \
     arm-20.1 \
     arm-20 \
     arm \
     default \
  ; then
  export ATDM_CONFIG_COMPILER=ARM-20.1_OPENMPI-4.0.5
else
  echo
  echo "***"
  echo "*** ERROR: A supported compiler was not selected for 'van1-tx2' (stria) env in buildname '${ATDM_CONFIG_BUILD_NAME}'"
  echo "***"
  echo "*** Supported compilers include:"
  echo "***"
  echo "****  arm-20.1-openmpi-4.0.5    (arm, arm-20, arm-20.1, default)"
  echo "***"
  return

fi
