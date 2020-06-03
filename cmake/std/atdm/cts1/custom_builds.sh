#
# Custom builds for 'cts1' env
#
# NOTE: This file gets sourced in atdm/utils/set_build_options.sh before the
# default grep logic is applied.
#

# Custom compiler selection logic

# TODO: Check with EMPIRE developers about default compiler/mpi ver
if atdm_match_any_buildname_keyword \
    intel-19.0.5-openmpi-4.0.1 \
    intel-19.0.5_openmpi-4.0.1 \
    intel-19.0.5 \
  ; then
  export ATDM_CONFIG_COMPILER=INTEL-19.0.5_OPENMPI-4.0.1
elif atdm_match_any_buildname_keyword \
      intel-19.0.4-openmpi-4.0.3 \
      intel-19.0.4_openmpi-4.0.3 \
      intel-19.0.4 \
      intel-19 \
      intel \
      default \
  ; then
  export ATDM_CONFIG_COMPILER=INTEL-19.0.4_OPENMPI-4.0.3
else
  echo
  echo "***"
  echo "*** ERROR: A supported compiler was not selected for 'cts1' env - $ATDM_CONFIG_BUILD_NAME"
  echo "***"
  echo "*** Supported compilers include:"
  echo "***"
  echo "****  intel-19.0.4_openmpi-4.0.3 (default)"
  echo "****  intel-19.0.5_openmpi-4.0.1"
  echo "***"
  return

fi
