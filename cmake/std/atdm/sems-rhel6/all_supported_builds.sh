# This script is sourced to return all of the supported builds

export ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX=Trilinos-atdm-

export ATDM_CONFIG_ALL_SUPPORTED_BUILDS=(
  sems-rhel6-gnu-debug-openmp
  sems-rhel6-gnu-debug-serial
  sems-rhel6-gnu-opt-openmp
  sems-rhel6-gnu-opt-serial
  sems-rhel6-intel-opt-openmp
  )
