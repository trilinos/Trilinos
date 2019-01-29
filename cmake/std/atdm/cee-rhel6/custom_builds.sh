#
# Custom builds for cee-rhel6 env
#
# NOTE: This file gets sourced in atdm/utils/set_build_options.sh before the
# default grep logic is applied.
#

# Custom compiler selection logic
if   [[ $ATDM_CONFIG_BUILD_NAME == *"clang-5.0.1-openmpi-1.10.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"clang-5.0.1_openmpi-1.10.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"clang-5.0.1"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"clang"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"default" ]] \
  ; then
  export ATDM_CONFIG_COMPILER=CLANG-5.0.1_OPENMPI-1.10.2

elif [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.2.0-openmpi-1.10.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.2.0_openmpi-1.10.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.2.0"* ]] \
  ; then
  export ATDM_CONFIG_COMPILER=GNU-7.2.0_OPENMPI-1.10.2

elif [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-4.9.3-openmpi-1.10.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-4.9.3_openmpi-1.10.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-4.9.3"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu"* ]] \
  ; then
  export ATDM_CONFIG_COMPILER=GNU-4.9.3_OPENMPI-1.10.2

elif [[ $ATDM_CONFIG_BUILD_NAME == *"intel-18.0.2-mpich2-3.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"intel-18.0.2_mpich2-3.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"intel-18.0.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"intel-18"* ]] \
  ; then
  export ATDM_CONFIG_COMPILER=INTEL-18.0.2_MPICH2-3.2

elif [[ $ATDM_CONFIG_BUILD_NAME == *"intel-17.0.1-intelmpi-5.1.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"intel-17.0.1_intelmpi-5.1.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"intel-17.0.1"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"intel-17"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"intel"* ]] \
  ; then
  export ATDM_CONFIG_COMPILER=INTEL-17.0.1_INTELMPI-5.1.2

else
  echo
  echo "***"
  echo "*** ERROR: A supported compiler was not selected for 'cee-rhel6' env"
  echo "***"
  echo "*** Supported compilers include:"
  echo "***"
  echo "****  clang-5.0.1-openmpi-1.10.2   (default)"
  echo "****  gnu-4.9.3-openmpi-1.10.2     (default gnu)"
  echo "****  gnu-7.2.0-openmpi-1.10.2"
  echo "****  intel-17.0.1-intelmpi-5.1.2  (default intel)"
  echo "****  intel-18.0.2-mpich2-3.2"
  echo "***"  
  return
fi
# ToDo: Add support for CUDA compilers above
