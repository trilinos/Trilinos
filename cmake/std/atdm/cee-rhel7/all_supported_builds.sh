# This script is sourced to return all of the supported builds

export ATDM_CONFIG_CTEST_S_BUILD_NAME_PREFIX=Trilinos-atdm-

HOSTNAME=$(hostname)

if [[ "${HOSTNAME}" == "ascicgpu"* ]] ; then
  export ATDM_CONFIG_ALL_SUPPORTED_BUILDS=(
    cee-rhel7_cuda-10.1.243_gnu-7.2.0_openmpi-4.0.3_shared_opt # SPARC Nightly build
    cee-rhel7_cuda-10.1.243_gnu-7.2.0_openmpi-4.0.3_shared_dbg # SPARC Nightly build
    ### SPARC MINI BUILDS ####
    cee-rhel7_mini_cuda-10.1.243_gnu-7.2.0_openmpi-4.0.3_static_opt
    cee-rhel7_mini_cuda-10.1.243_gnu-7.2.0_openmpi-4.0.3_static_dbg
    cee-rhel7_mini-no-mpi_cuda-10.1.243_gnu-7.2.0_static_opt
    )
else
  export ATDM_CONFIG_ALL_SUPPORTED_BUILDS=(
    #cee-rhel7_clang-9.0.1_openmpi-4.0.3_serial_static_dbg     # SPARC has installs with this build
    cee-rhel7_clang-9.0.1_openmpi-4.0.3_serial_static_opt      # SPARC CI build
    cee-rhel7_gnu-7.2.0_openmpi-4.0.3_serial_shared_opt        # SPARC CI build
    cee-rhel7_intel-19.0.3_mpich2-3.2_openmp_static_opt        # SPARC CI build
    cee-rhel7_intel-19.0.3_intelmpi-2018.4_serial_static_opt   # SPARC Nightly bulid
    cee-rhel7_mini_clang-9.0.1_openmpi-4.0.3_serial_static_opt
    #### SPARC MINI BUILDS ####
    cee-rhel7_mini_gnu-7.2.0_openmpi-4.0.3_serial_static_opt
    cee-rhel7_mini_clang-9.0.1_openmpi-4.0.3_serial_static_opt
    cee-rhel7_mini-no-mpi_gnu-7.2.0_serial_static_opt
    cee-rhel7_mini-no-mpi_clang-9.0.1_serial_static_opt
    cee-rhel7_mini_intel-19.0.3_mpich2-3.2_static_opt
    cee-rhel7_mini-no-mpi_intel-19.0.3_static_opt
    )
fi

# NOTE: Above, we have commented out the 'dbg' build because it was running
# the test suite very slow and had many timeouts (see ATDV-322)

# NOTE: You have to run the CEE CUDA builds on an 'ascicgpu' machine or they
# generate binaries that don't work for the 'ascicgpu' machines.  At the same
# time, the 'ascicgpu' machines are a coveted resource so we don't want to
# waste an 'ascicgpu' machine with builds that don't even require running on
# the GPU.  So to run the non CUDA-builds, get on a CEE machine that does not
# have a GPU and to run the CUDA builds, get on an 'ascicgpu' machine.
