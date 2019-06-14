# This script is sourced to return all of the supported builds

export ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX=Trilinos-atdm-

export ATDM_CONFIG_ALL_SUPPORTED_BUILDS=(
  spack-rhel-gnu-7.2.0-openmp-debug
  spack-rhel-gnu-7.2.0-openmp-release
  spack-rhel-gnu-7.2.0-openmp-release-debug
  spack-rhel-gnu-7.2.0-serial-debug
  spack-rhel-gnu-7.2.0-serial-release
  spack-rhel-gnu-7.2.0-serial-release-debug
  )
