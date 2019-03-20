# This script is sourced to return all of the supported builds

export ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX=Trilinos-atdm-

export ATDM_CONFIG_ALL_SUPPORTED_BUILDS=(
  spack-rhel6-gnu-debug-openmp
  spack-rhel6-gnu-debug-serial
  spack-rhel6-gnu-opt-openmp
  spack-rhel6-gnu-opt-serial
  )
