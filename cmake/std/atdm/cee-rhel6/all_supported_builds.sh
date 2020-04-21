# This script is sourced to return all of the supported builds

export ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX=Trilinos-atdm-

export ATDM_CONFIG_ALL_SUPPORTED_BUILDS=(
  #cee-rhel6_clang-9.0.1_openmpi-4.0.2_serial_static_dbg    # SPARC has installs with this build
  cee-rhel6_clang-9.0.1_openmpi-4.0.2_serial_static_opt    # SPARC CI build
  cee-rhel6_gnu-7.2.0_openmpi-4.0.2_serial_shared_opt      # SPARC CI build
  cee-rhel6_intel-18.0.2_mpich2-3.2_openmp_static_opt      # SPARC CI build
  cee-rhel6_intel-19.0.3_intelmpi-2018.4_serial_static_opt # SPARC Nightly bulid
  )

# NOTE: Above, we have commented out the 'dbg' build because it was running
# the test suite very slow and had many timeouts (see ATDV-322)
