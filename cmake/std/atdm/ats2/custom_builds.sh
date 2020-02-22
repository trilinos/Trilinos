#
# Custom builds for 'ats2' env
#
# NOTE: This file gets sourced in atdm/utils/set_build_options.sh before the
# default grep logic is applied.
#

# Custom compiler selection logic

if   [[ $ATDM_CONFIG_BUILD_NAME == *"cuda-10.1.243-xl-2019.08.20-spmpi-2019.06.24"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"cuda-10.1.243_xl-2019.08.20_spmpi-2019.06.24"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"cuda-10.1.243-xl-2019.08.20"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"cuda-10.1.243_xl-2019.08.20"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"cuda-10.1.243-xl-2019"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"cuda-10.1.243_xl-2019"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"cuda-xl"* ]] \
  ; then
  export ATDM_CONFIG_COMPILER=CUDA-10.1.243_XL-2019.08.20_SPMPI-2019.06.24

elif [[ $ATDM_CONFIG_BUILD_NAME == *"xl-2019.08.20-spmpi-2019.06.24"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"xl-2019.08.20_spmpi-2019.06.24"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"xl-2019.08.20"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"xl-2019"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"xl"* ]] \
  ; then
  export ATDM_CONFIG_COMPILER=XL-2019.08.20_SPMPI-2019.06.24

elif [[ $ATDM_CONFIG_BUILD_NAME == *"cuda-10.1.243-gnu-7.3.1-spmpi-2019.06.24"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"cuda-10.1.243_gnu-7.3.1_spmpi-2019.06.24"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"cuda-10.1.243-gnu-7.3.1"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"cuda-10.1.243_gnu-7.3.1"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"cuda-10.1.243-gnu-7"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"cuda-10.1.243_gnu-7"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"cuda-10.1.243"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"cuda-10"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"cuda-gnu"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"cuda"* ]] \
  ; then
  export ATDM_CONFIG_COMPILER=CUDA-10.1.243_GNU-7.3.1_SPMPI-2019.06.24
  # NOTE: Default 'cuda' must be last cuda listed!

elif [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.3.1-spmpi-2019.06.24"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.3.1_spmpi-2019.06.24"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.3.1"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"default" ]] \
  ; then
  export ATDM_CONFIG_COMPILER=GNU-7.3.1_SPMPI-2019.06.24
  # NOTE: Defaut 'gnu' must be last 'gnu' listed!

else
  echo
  echo "***"
  echo "*** ERROR: A supported compiler was not selected for 'ats2' env - $ATDM_CONFIG_BUILD_NAME"
  echo "***"
  echo "*** Supported compilers include:"
  echo "***"
  echo "****  gnu-7.3.1_spmpi-2019.06.24                  (default, default gnu)"
  echo "****  cuda-10.1.243_gnu-7.3.1_spmpi-2019.06.24    (default cuda)"
  echo "****  xl-2019.08.20_spmpi-2019.06.24              (disabled)"
  echo "****  cuda-10.1.243-gnu-7.3.1-spmpi-2019.06.24    (disabled)"
  echo "***"
  return

fi
