#
# Custom builds for cee-rhel6 env
#
# NOTE: This file gets sourced in atdm/utils/set_build_options.sh before the
# default grep logic is applied.
#

# Custom compiler selection logic

if   [[ $ATDM_CONFIG_BUILD_NAME == *"clang-5.0.1-openmpi-1.10.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"clang-5.0.1_openmpi-1.10.2"* ]] \
  ; then
  export ATDM_CONFIG_COMPILER=CLANG-5.0.1_OPENMPI-1.10.2

elif [[ $ATDM_CONFIG_BUILD_NAME == *"clang-5.0.1-openmpi-4.0.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"clang-5.0.1_openmpi-4.0.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"clang-5.0.1"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"clang-5"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"clang"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"default" ]] \
  ; then
  export ATDM_CONFIG_COMPILER=CLANG-5.0.1_OPENMPI-4.0.2
  # Must list the default clang build last for correct matching of of defaults

elif [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.2.0-openmpi-1.10.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.2.0_openmpi-1.10.2"* ]] \
  ; then
  export ATDM_CONFIG_COMPILER=GNU-7.2.0_OPENMPI-1.10.2

elif [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.2.0-openmpi-4.0.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.2.0_openmpi-4.0.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7.2.0"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu-7"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"gnu"* ]] \
  ; then
  export ATDM_CONFIG_COMPILER=GNU-7.2.0_OPENMPI-4.0.2
  # List default "gnu"* build last for correct matching of defaults

elif [[ $ATDM_CONFIG_BUILD_NAME == *"intel-18.0.2-mpich2-3.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"intel-18.0.2_mpich2-3.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"intel-18.0.2"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"intel-18"* ]] \
  ; then
  export ATDM_CONFIG_COMPILER=INTEL-18.0.2_MPICH2-3.2

elif [[ $ATDM_CONFIG_BUILD_NAME == *"intel-19.0.3-intelmpi-2018.4"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"intel-19.0.3_intelmpi-2018.4"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"intel-19.0.3"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"intel-19"* ]] \
  || [[ $ATDM_CONFIG_BUILD_NAME == *"intel"* ]] \
  ; then
  export ATDM_CONFIG_COMPILER=INTEL-19.0.3_INTELMPI-2018.4
  # List default "intel"* build last for correct matching!

else
  echo
  echo "***"
  echo "*** ERROR: A supported compiler was not selected for 'cee-rhel6' env"
  echo "***"
  echo "*** Supported compilers include:"
  echo "***"
  echo "****  clang-5.0.1-openmpi-1.10.2     (DEPRECATED)"
  echo "****  clang-5.0.1-openmpi-4.0.2      (default, default clang)"
  echo "****  gnu-7.2.0-openmpi-1.10.2       (DEPRECATED)"
  echo "****  gnu-7.2.0-openmpi-4.0.2        (default gnu)"
  echo "****  intel-18.0.2-mpich2-3.2        (DEPRECATED)"
  echo "****  intel-19.0.3-intelmpi-2018.4   (default intel)"
  echo "***"  
  return

fi
# ToDo: Add support for CUDA compilers above
