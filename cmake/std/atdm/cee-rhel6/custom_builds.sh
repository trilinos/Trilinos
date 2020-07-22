#
# Custom builds for cee-rhel6 env
#
# NOTE: This file gets sourced in atdm/utils/set_build_options.sh before the
# default grep logic is applied.
#

# Custom compiler selection logic

if atdm_match_any_buildname_keyword \
    clang-9.0.1-openmpi-4.0.3 \
    clang-9.0.1_openmpi-4.0.3 \
    clang-9.0.1 \
    clang-9 \
    clang \
    default \
  ; then
  export ATDM_CONFIG_COMPILER=CLANG-9.0.1_OPENMPI-4.0.3
  # Must list the default clang build last of all the 'clang' builds for
  # correct matching of of defaults

elif atdm_match_any_buildname_keyword \
    gnu-7.2.0-openmpi-4.0.3 \
    gnu-7.2.0_openmpi-4.0.3 \
    gnu-7.2.0 \
    gnu-7 \
    gnu \
  ; then
  export ATDM_CONFIG_COMPILER=GNU-7.2.0_OPENMPI-4.0.3
  # List default "gnu"* build last of all the 'gnu' builds for correct
  # matching of defaults

elif atdm_match_any_buildname_keyword \
    intel-18.0.2-mpich2-3.2 \
    intel-18.0.2_mpich2-3.2 \
    intel-18.0.2 \
    intel-18 \
  ; then
  export ATDM_CONFIG_COMPILER=INTEL-18.0.2_MPICH2-3.2

elif atdm_match_any_buildname_keyword \
    intel-19.0.3-intelmpi-2018.4 \
    intel-19.0.3_intelmpi-2018.4 \
    intel-19.0.3 \
    intel-19 \
    intel \
  ; then
  export ATDM_CONFIG_COMPILER=INTEL-19.0.3_INTELMPI-2018.4
  # List default intel build last of all the 'intel' builds for correct
  # matching!

elif atdm_match_any_buildname_keyword \
    cuda-10.1.243_gcc-7.2.0-openmpi-4.0.3 \
    cuda-10.1.243_gcc-7.2.0_openmpi-4.0.3 \
    cuda-10.1.243_gcc-7.2.0 \
    cuda-10.1.243 \
    cuda-10 \
    cuda \
  ; then
  export ATDM_CONFIG_COMPILER=CUDA-10.1.243_GCC-7.2.0_OPENMPI-4.0.3

else
  echo
  echo "***"
  echo "*** ERROR: A supported compiler was not selected for 'cee-rhel6' env in buildname '${ATDM_CONFIG_BUILD_NAME}'"
  echo "***"
  echo "*** Supported compilers include:"
  echo "***"
  echo "****  clang-9.0.1-openmpi-4.0.3               (default, default clang)"
  echo "****  gnu-7.2.0-openmpi-4.0.3                 (default gnu)"
  echo "****  intel-18.0.2-mpich2-3.2"
  echo "****  intel-19.0.3-intelmpi-2018.4            (default intel)"
  echo "****  cuda-10.1.243_gcc-7.2.0-openmpi-4.0.3   (default cuda)"
  echo "***"  
  return

fi
# ToDo: Add support for CUDA compilers above
