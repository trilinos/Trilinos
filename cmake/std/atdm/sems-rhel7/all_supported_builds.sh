# This script is sourced to return all of the supported builds

export ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX=Trilinos-atdm-

export ATDM_CONFIG_ALL_SUPPORTED_BUILDS=(
  sems-rhel7-cuda-9.2-Volta70-complex-shared-release-debug
  sems-rhel7-gnu-7.2.0-openmp-complex-shared-release-debug
  sems-rhel7-intel-18.0.5-openmp-complex-shared-release-debug
  sems-rhel7-intel-18.0.5-openmp-shared-debug
  sems-rhel7-intel-18.0.5-openmp-shared-release-debug
  )
