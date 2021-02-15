# This file contains the options needed to both run the pull request testing
# for Trilinos for the Linux Clang 10.0.0 pull request testing builds, and to reproduce
# the errors reported by those builds. Prior to using this this file, the
# appropriate set of SEMS modules must be loaded and accessible through the
# SEMS NFS mount. (See the sems/PullRequest*TestingEnv.sh files.)

# Usage: cmake -C PullRequestLinuxClang10.0.0TestingSettings.cmake

# Misc options typically added by CI testing mode in TriBITS

# Use the below option only when submitting to the dashboard
#set (CTEST_USE_LAUNCHERS ON CACHE BOOL "Set by default for PR testing")

#set (TPL_ENABLE_Netcdf OFF CACHE BOOL "Turn off for Clang")

#set (TPL_Netcdf_LIBRARIES "-L$ENV{SEMS_NETCDF_ROOT}/lib;-L$ENV{SEMS_HDF5_ROOT}/lib;$ENV{SEMS_NETCDF_ROOT}/lib/libnetcdf.a;$ENV{SEMS_NETCDF_ROOT}/lib/libpnetcdf.a;$ENV{SEMS_HDF5_ROOT}/lib/libhdf5_hl.a;$ENV{SEMS_HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl;-lcurl" CACHE STRING "Set by default for CUDA PR testing")

set (MPI_EXEC_PRE_NUMPROCS_FLAGS "--bind-to;none" CACHE STRING "Set by default for PR testing")
# NOTE: The above is a workaround for the problem of having threads on MPI
# ranks bind to the same cores (see #2422).

# Disable just one Teko sub-unit test that fails with openmpi 1.10 (#2712)
set (Teko_DISABLE_LSCSTABALIZED_TPETRA_ALPAH_INV_D ON CACHE BOOL "Temporarily disabled in PR testing")

# Disable three ShyLu_DD tests - see #2691
set (ShyLU_DDFROSch_test_frosch_laplacian_epetra_2d_gdsw_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ShyLU_DDFROSch_test_frosch_laplacian_epetra_2d_rgdsw_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ShyLU_DDFROSch_test_frosch_interfacesets_2D_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")

include("${CMAKE_CURRENT_LIST_DIR}/PullRequestLinuxCommonTestingSettings.cmake")

#Disable for clang
set(FEI_elemDOF_Aztec_MPI_2_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set(FEI_lagrange_20quad_old_MPI_2_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set(FEI_lagrange_20quad_old_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set(FEI_multifield_vbr_az_MPI_2_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set(FEI_multifield_vbr_az_MPI_3_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set(ROL_example_PinT_parabolic-control_example_01_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set(Rythmos_StepperBuilder_UnitTest_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
