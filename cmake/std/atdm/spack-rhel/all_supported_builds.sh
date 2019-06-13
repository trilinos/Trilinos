# This script is sourced to return all of the supported builds

export ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX=Trilinos-atdm-

export ATDM_CONFIG_ALL_SUPPORTED_BUILDS=(
  spack-rhel-gnu-debug-openmp
  spack-rhel-gnu-debug-serial
  spack-rhel-gnu-opt-openmp
  spack-rhel-gnu-opt-serial
  )
