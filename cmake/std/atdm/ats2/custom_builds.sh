#
# Custom builds for 'ats2' env
#
# NOTE: This file gets sourced in atdm/utils/set_build_options.sh before the
# default grep logic is applied.
#

# Custom compiler selection logic
if atdm_match_any_buildname_keyword \
  cuda-10.1.243-xl-2020.03.18-spmpi-rolling \
  cuda-10.1.243-xl-2020.03.18_spmpi-rolling \
  cuda-10.1.243_xl-2020.03.18-spmpi-rolling \
  cuda-10.1.243_xl-2020.03.18_spmpi-rolling \
  cuda-10.1.243-xl-2020.03.18 \
  cuda-10.1.243_xl-2020.03.18 \
  cuda-10.1.243-xl-2020 \
  cuda-10.1.243_xl-2020 \
  cuda-xl \
  cuda_xl \
  ; then
  export ATDM_CONFIG_COMPILER=CUDA-10.1.243_XL-2020.03.18_SPMPI-ROLLING

elif atdm_match_any_buildname_keyword \
  xl-2020.03.18-spmpi-rolling \
  xl-2020.03.18_spmpi-rolling \
  xl-2020.03.18 \
  xl-2020 \
  xl \
  ; then
  export ATDM_CONFIG_COMPILER=XL-2020.03.18_SPMPI-ROLLING

elif atdm_match_any_buildname_keyword \
  cuda-10.1.243-gnu-7.3.1-spmpi-rolling \
  cuda-10.1.243_gnu-7.3.1_spmpi-rolling \
  cuda-10.1.243-gnu-7.3.1 \
  cuda-10.1.243_gnu-7.3.1 \
  cuda-10.1.243-gnu-7 \
  cuda-10.1.243_gnu-7 \
  cuda-10.1.243 \
  cuda-10 \
  cuda-gnu \
  cuda \
  ; then
  export ATDM_CONFIG_COMPILER=CUDA-10.1.243_GNU-7.3.1_SPMPI-ROLLING
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
  echo "****  cuda-10.1.243_gnu-7.3.1_spmpi-rolling        (default cuda)"
  echo "****  xl-2020.03.18_spmpi-rolling                  (default xl)"
  echo "****  cuda-10.1.243-xl-2020.03.18_spmpi-rolling    (default cuda-xl)"
  echo "***"
  return

fi
