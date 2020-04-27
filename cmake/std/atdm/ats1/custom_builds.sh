#
# Custom builds for 'ats1' env
#
# NOTE: This file gets sourced in atdm/utils/set_build_options.sh before the
# default grep logic is applied.
#

# Custom compiler selection logic

if atdm_match_any_buildname_keyword \
     intel-19.0.4-mpich-7.7.6 \
     intel-19.0.4_mpich-7.7.6 \
     intel-19.0.4 \
     intel-19 \
  ; then
  export ATDM_CONFIG_COMPILER=INTEL-19.0.4_MPICH-7.7.6
# Note: matching "intel" below is dangerous. The user could easily make a typo
# and wind up with the intel-18 env. Leaving "intel" as a match below due to
# customer request.
elif atdm_match_any_buildname_keyword \
       intel-18.0.5-mpich-7.7.6 \
       intel-18.0.5_mpich-7.7.6 \
       intel-18.0.5 \
       intel-18 \
       default \
       intel \
  ; then
  export ATDM_CONFIG_COMPILER=INTEL-18.0.5_MPICH-7.7.6
else
  echo
  echo "***"
  echo "*** ERROR: A supported compiler was not selected for 'ats1' env - $ATDM_CONFIG_BUILD_NAME"
  echo "***"
  echo "*** Supported compilers include:"
  echo "***"
  echo "****  intel-19.0.4-mpich-7.7.6                  (intel-19)"
  echo "****  intel-18.0.5-mpich-7.7.6                  (intel-18, intel, default)"
  echo "***"
  return

fi
