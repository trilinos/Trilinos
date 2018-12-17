# This script is sourced to return all of the supported builds

export ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX=Trilinos-atdm-

export ATDM_CONFIG_ALL_SUPPORTED_BUILDS=(
  sems-rhel7-clang-3.9.0-openmp-release-debug
  sems-rhel7-clang-3.9.0-serial-release-debug
  sems-rhel7-gnu-7.2.0-openmp-release-debug
  sems-rhel7-gnu-7.2.0-serial-release-debug
  #sems-rhel7-intel-17.0.1-openmp-release-debug
  #sems-rhel7-intel-17.0.1-serial-release-debug
  sems-rhel7-cuda-9.2-Pascal60-opt
  sems-rhel7-cuda-9.2-Pascal60-release-debug
  )

# NOTE: As of 12/8/2018, the intel-17.0.1 builds are too slow to practically
# run (something wrong with the Intel compiler setup).
