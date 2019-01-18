# This file contains the options needed to both run the pull request testing
# for Trilinos for the CUDA 9.2 pull request testing builds, and to reproduce
# the errors reported by those builds. Prior to using this this file, the
# appropriate set of modules must be loaded and path must be augmented.
# (See the sems/PullRequestCuda9.2TestingEnv.sh files.)

# Usage: cmake -C PullRequestLinuxCUDA9.2TestingSettings.cmake

# Misc options typically added by CI testing mode in TriBITS

# Use the below option only when submitting to the dashboard
set (CTEST_USE_LAUNCHERS ON CACHE BOOL "Set by default for PR testing")

# Options necessary for CUDA build
set (TPL_ENABLE_MPI ON CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_CUDA ON CACHE BOOL "Set by default for CUDA PR testing")
set (Kokkos_ENABLE_Cuda ON CACHE BOOL "Set by default for CUDA PR testing")
set (Kokkos_ENABLE_Cuda_UVM ON CACHE BOOL "Set by default for CUDA PR testing")
set (KOKKOS_ARCH Power8 CACHE STRING "Set by default for CUDA PR testing")

# TPL settings specific to CUDA build
set (TPL_BLAS_LIBRARIES "-L${BLAS_ROOT}/lib -lblas -lgfortran -lgomp -lm" CACHE STRING "Set by default for CUDA PR testing")
set (TPL_LAPACK_LIBRARIES "-L${LAPACK_ROOT}/lib -llapack -lgfortran -lgomp" CACHE STRING "Set by default for CUDA PR testing")
set (BUILD_SHARED_LIBS OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Tpetra_INST_SERIAL ON CACHE BOOL "Set by default for CUDA PR testing")
set (Trilinos_ENABLE_SECONDARY_TESTED_CODE OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_Scotch OFF CACHE BOOL "Set by default for CUDA PR testing")
# Parmetis is available on ride and could be enabled for the CUDA PR build
set (TPL_ENABLE_ParMETIS OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_Netcdf_LIBRARIES "-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;-L${HDF5_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.a;${BOOST_ROOT}/lib/libboost_system.a;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl" CACHE STRING "Set by default for CUDA PR testing")
# SuperLU is available on ride and could be enabled for the CUDA PR build
set (TPL_ENABLE_SuperLU OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_BoostLib OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Trilinos_ENABLE_Moertel OFF CACHE BOOL "Disable for CUDA PR testing")
set (Trilinos_ENABLE_Komplex OFF CACHE BOOL "Disable for CUDA PR testing")

# Temporary options to clean up build
set (Trilinos_ENABLE_SEACAS OFF CACHE BOOL "Temporary disable for CUDA PR testing")
set (Trilinos_ENABLE_DEBUG OFF CACHE BOOL "Temporary disable for CUDA PR testing")
set (Trilinos_ENABLE_DEBUG_SYMBOLS OFF CACHE BOOL "Temporary disable for CUDA PR testing")
set (STK_ENABLE_TESTS OFF CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_interfacepartitionofunity_DIM2_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_interfacesets_DIM2_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_localpartitionofunitybasis_EPETRA_MPI_1_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_GDSWP_DIM2_DPN1_ORD0_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_GDSWP_DIM2_DPN1_ORD1_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_GDSWP_DIM2_DPN2_ORD0_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_GDSWP_DIM2_DPN2_ORD1_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_RGDSWP_DIM2_DPN1_ORD0_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_RGDSWP_DIM2_DPN1_ORD1_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_RGDSWP_DIM2_DPN2_ORD0_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_RGDSWP_DIM2_DPN2_ORD1_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM2_DPN1_ORD0_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM2_DPN1_ORD1_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM2_DPN2_ORD0_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM2_DPN2_ORD1_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM2_DPN1_ORD0_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM2_DPN1_ORD1_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM2_DPN2_ORD0_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM2_DPN2_ORD1_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_GDSW_DIM2_DPN1_ORD0_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_GDSW_DIM2_DPN1_ORD1_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_GDSW_DIM2_DPN2_ORD0_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_GDSW_DIM2_DPN2_ORD1_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_IPOUHarmonic_DIM2_DPN1_ORD0_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_IPOUHarmonic_DIM2_DPN1_ORD1_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_IPOUHarmonic_DIM2_DPN2_ORD0_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_IPOUHarmonic_DIM2_DPN2_ORD1_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_RGDSW_DIM2_DPN1_ORD0_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_RGDSW_DIM2_DPN1_ORD1_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_RGDSW_DIM2_DPN2_ORD0_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_RGDSW_DIM2_DPN2_ORD1_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_stokes_hdf5_TLBP_GDSW_O0_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ShyLU_DDFROSch_test_thyra_xpetra_stokes_hdf5_TLBP_GDSW_O1_EPETRA_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ROL_example_PDE-OPT_0ld_adv-diff-react_example_01_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ROL_example_PDE-OPT_0ld_adv-diff-react_example_02_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ROL_example_PDE-OPT_0ld_poisson_example_01_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ROL_example_PDE-OPT_0ld_stefan-boltzmann_example_03_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ROL_example_PDE-OPT_navier-stokes_example_01_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ROL_example_PDE-OPT_navier-stokes_example_02_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ROL_example_PDE-OPT_nonlinear-elliptic_example_01_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ROL_example_PDE-OPT_nonlinear-elliptic_example_02_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ROL_example_PDE-OPT_obstacle_example_01_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ROL_example_PDE-OPT_topo-opt_poisson_example_01_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ROL_test_elementwise_TpetraMultiVector_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (TpetraCore_Core_initialize_where_tpetra_initializes_kokkos_MPI_1_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (TpetraCore_Core_initialize_where_tpetra_initializes_mpi_and_user_initializes_kokkos_MPI_2_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (TpetraCore_Core_initialize_where_user_initializes_kokkos_MPI_1_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (TpetraCore_Core_initialize_where_user_initializes_mpi_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (TpetraCore_Core_ScopeGuard_where_tpetra_initializes_kokkos_MPI_1_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (TpetraCore_Core_ScopeGuard_where_tpetra_initializes_mpi_and_user_initializes_kokkos_MPI_2_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (TpetraCore_Core_ScopeGuard_where_user_initializes_kokkos_MPI_1_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (TpetraCore_Core_ScopeGuard_where_user_initializes_mpi_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (TrilinosCouplings_Example_Maxwell_MueLu_MPI_1_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (TrilinosCouplings_Example_Maxwell_MueLu_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")

include("${CMAKE_CURRENT_LIST_DIR}/PullRequestLinuxCommonTestingSettings.cmake")

