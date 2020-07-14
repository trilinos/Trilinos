#
# Custom builds for van1-tx2 env
#

# Try matching against arm-20.1 before arm-20 or arm
if atdm_match_any_buildname_keyword \
     arm-20.1-openmpi-4.0.3 \
     arm-20.1_openmpi-4.0.3 \
     arm-20.1 \
  ; then
  export ATDM_CONFIG_COMPILER=ARM-20.1_OPENMPI-4.0.3
elif atdm_match_any_buildname_keyword \
     arm-20.0-openmpi-4.0.2 \
     arm-20.0_openmpi-4.0.2 \
     arm-20.0 \
     arm-20 \
     arm \
     default \
  ; then
  export ATDM_CONFIG_COMPILER=ARM-20.0_OPENMPI-4.0.2
else
  echo
  echo "***"
  echo "*** ERROR: A supported compiler was not selected for 'van1-tx2' (stria) env in buildname '${ATDM_CONFIG_BUILD_NAME}'"
  echo "***"
  echo "*** Supported compilers include:"
  echo "***"
  echo "****  arm-20.0-openmpi-4.0.2    (arm-20.0, default)"
  echo "****  arm-20.1-openmpi-4.0.3    (arm-20.1)"
  echo "***"
  return

fi
