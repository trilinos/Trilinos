#
# Custom builds for 'ats2' env
#
# NOTE: This file gets sourced in atdm/utils/set_build_options.sh before the
# default grep logic is applied.
#

# Custom compiler selection logic
if atdm_match_any_buildname_keyword \
  cuda-11.2.152-gcc-8.3.1-spmpi-rolling \
  cuda-11.2.152_gcc-8.3.1_spmpi-rolling \
  cuda-11.2.152-gcc-8.3.1 \
  cuda-11.2.152_gcc-8.3.1 \
  cuda-11.2.152-gcc-8 \
  cuda-11.2.152_gcc-8 \
  cuda-11.2.152 \
  cuda-11 \
  cuda-gnu \
  cuda \
  ; then
  export ATDM_CONFIG_COMPILER=CUDA-11.2.152_GCC-8.3.1_SPMPI-ROLLING
  # NOTE: Default 'cuda' must be last cuda listed!

elif atdm_match_any_buildname_keyword \
  gnu-7.3.1-spmpi-rolling \
  gnu-7.3.1_spmpi-rolling \
  gnu-7.3.1 \
  gnu-7 \
  gnu \
  default \
  ; then
  export ATDM_CONFIG_COMPILER=GNU-7.3.1_SPMPI-ROLLING
  # NOTE: Defaut 'gnu' must be last 'gnu' listed!

else
  echo
  echo "***"
  echo "*** ERROR: A supported compiler was not selected for 'ats2' env - $ATDM_CONFIG_BUILD_NAME"
  echo "***"
  echo "*** Supported compilers include:"
  echo "***"
  echo "****  gnu-7.3.1_spmpi-rolling                      (default, default gnu)"
  echo "****  cuda-11.2.152_gcc-8.3.1_spmpi-rolling        (default cuda)"
  echo "***"
  return

fi
