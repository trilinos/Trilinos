#
# Custom builds for 'cts1empire' env
#
# NOTE: This file gets sourced in atdm/utils/set_build_options.sh before the
# default grep logic is applied.
#

# Custom compiler selection logic

if atdm_match_any_buildname_keyword \
    intel-18.0.2-openmpi-4.0.1 \
    intel-18.0.2_openmpi-4.0.1 \
    intel-18.0.2 \
    intel \
    default \
  ; then
  export ATDM_CONFIG_COMPILER=INTEL-18.0.2_OPENMPI-4.0.1
else
  echo
  echo "***"
  echo "*** ERROR: A supported compiler was not selected for 'cts1empire' env - $ATDM_CONFIG_BUILD_NAME"
  echo "***"
  echo "*** Supported compilers include:"
  echo "***"
  echo "****  intel-18.0.2_openmpi-4.0.1   (default)"
  echo "***"
  return

fi
