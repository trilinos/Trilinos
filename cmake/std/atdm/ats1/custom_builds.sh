#
# Custom builds for 'ats1' env
#
# NOTE: This file gets sourced in atdm/utils/set_build_options.sh before the
# default grep logic is applied.
#

# Custom compiler selection logic

if atdm_match_any_buildname_keyword \
     intel-19.0.4-mpich-7.7.15 \
     intel-19.0.4_mpich-7.7.15 \
     intel-19.0.4 \
     intel-19 \
     intel \
     default \
  ; then
  export ATDM_CONFIG_COMPILER=INTEL-19.0.4_MPICH-7.7.15
else
  echo
  echo "***"
  echo "*** ERROR: A supported compiler was not selected for 'ats1' env - $ATDM_CONFIG_BUILD_NAME"
  echo "***"
  echo "*** Supported compilers include:"
  echo "***"
  echo "****  intel-19.0.4-mpich-7.7.15                  (intel-19, intel, default)"
  echo "***"
  return

fi
