# This file contains the options needed to both run the pull request testing
# for Trilinos for the Linux Clang 11.0.1 pull request testing builds, and to reproduce
# the errors reported by those builds. Prior to using this this file, the
# appropriate set of SEMS modules must be loaded and accessible through the
# SEMS NFS mount. (See the sems/PullRequest*TestingEnv.sh files.)

# Usage: cmake -C PullRequestLinuxClang11.0.1TestingSettings.cmake

# Misc options typically added by CI testing mode in TriBITS

# Use the below option only when submitting to the dashboard
#set (CTEST_USE_LAUNCHERS ON CACHE BOOL "Set by default for PR testing")

#set (TPL_ENABLE_Netcdf OFF CACHE BOOL "Turn off for Clang")

#set (TPL_Netcdf_LIBRARIES "-L$ENV{SEMS_NETCDF_ROOT}/lib;-L$ENV{SEMS_HDF5_ROOT}/lib;$ENV{SEMS_NETCDF_ROOT}/lib/libnetcdf.a;$ENV{SEMS_NETCDF_ROOT}/lib/libpnetcdf.a;$ENV{SEMS_HDF5_ROOT}/lib/libhdf5_hl.a;$ENV{SEMS_HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl;-lcurl" CACHE STRING "Set by default for CUDA PR testing")

set (MPI_EXEC_PRE_NUMPROCS_FLAGS "--bind-to;none" CACHE STRING "Set by default for PR testing")
# NOTE: The above is a workaround for the problem of having threads on MPI
# ranks bind to the same cores (see #2422).

include("${CMAKE_CURRENT_LIST_DIR}/PullRequestLinuxCommonTestingSettings.cmake")

#Disable failing test in new Clang 11.0.1 build
set(Zoltan_ch_simple_parmetis_parallel_DISABLE ON CACHE BOOL "Disabled in PR testing")
