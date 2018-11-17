# This script is sourced to return all of the supported builds

export ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX=Trilinos-atdm-

export ATDM_CONFIG_ALL_SUPPORTED_BUILDS=(
  cee-rhel6-clang-5.0.1-openmpi-1.10.2-serial-static-opt
  cee-rhel6-gnu-4.9.3-openmpi-1.10.2-serial-static-opt
  cee-rhel6-gnu-7.2.0-openmpi-1.10.2-serial-static-opt
  cee-rhel6-intel-17.0.1-intelmpi-5.1.2-serial-static-opt
  cee-rhel6-intel-18.0.2-mpich2-3.2-serial-static-opt
  )
