#
# Custom builds for 'ats2' env
#
# NOTE: This file gets sourced in atdm/utils/set_build_options.sh before the
# default grep logic is applied.
#

# Custom compiler selection logic

# TODO: Check with SPARC developers about default compiler/mpi versions
if [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.3.1-spmpi-2019.06.24"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.3.1_spmpi-2019.06.24"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.3.1"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7"* ]] \
  ; then
  export ATDM_CONFIG_COMPILER=GNU-7.3.1_SPMPI-2019.06.24
else
  echo
  echo "***"
  echo "*** ERROR: A supported compiler was not selected for 'ats2' env"
  echo "***"
  echo "*** Supported compilers include:"
  echo "***"
  echo "****  gnu-7.3.1_spmpi-2019.06.24   (default)"
  echo "***"
  return

fi
