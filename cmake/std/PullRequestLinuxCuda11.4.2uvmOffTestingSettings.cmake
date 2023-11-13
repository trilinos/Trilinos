# This file contains the options needed to both run the pull request testing
# for Trilinos for the CUDA 10.1.105 UVM pull request testing builds, and to reproduce
# the errors reported by those builds. Prior to using this this file, the
# appropriate set of modules must be loaded and path must be augmented.
# (See the sems/PullRequestCuda10.1.105uvmOffTestingEnv.sh files.)

# Usage: cmake -C PullRequestLinuxCUDA10.1.105uvmOffTestingSettings.cmake

set (CMAKE_CXX_STANDARD "14" CACHE STRING "Set C++ standard to C++14")

# Misc options typically added by CI testing mode in TriBITS

# Use the below option only when submitting to the dashboard
set (CTEST_USE_LAUNCHERS ON CACHE BOOL "Set by default for PR testing")

# Options necessary for CUDA build
set (TPL_ENABLE_MPI ON CACHE BOOL "Set by default for CUDA PR testing")
set (Kokkos_ENABLE_CUDA ON CACHE BOOL "Set by default for CUDA PR testing")
set (KOKKOS_ARCH VOLTA70 CACHE STRING "Set by default for CUDA PR testing")
set (Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Sacado_ENABLE_HIERARCHICAL_DFAD ON CACHE BOOL "Set by default for CUDA PR testing")
set (Kokkos_ENABLE_CXX11_DISPATCH_LAMBDA ON CACHE BOOL "Set by default for CUDA PR testing")
set (Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "Set by default for CUDA PR testing")
set (Phalanx_KOKKOS_DEVICE_TYPE CUDA CACHE STRING "Set by default for CUDA PR testing")
set (MPI_EXEC_POST_NUMPROCS_FLAGS "-map-by;socket:PE=4" CACHE STRING "Set by default for CUDA PR testing")

# Options set to match the ATDM build
set (Trilinos_ENABLE_DEBUG OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Trilinos_ENABLE_DEBUG_SYMBOLS OFF CACHE BOOL "Set by default for CUDA PR testing")
set (BUILD_SHARED_LIBS OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Tpetra_INST_SERIAL ON CACHE BOOL "Set by default for CUDA PR testing")
set (Trilinos_ENABLE_SECONDARY_TESTED_CODE OFF CACHE BOOL "Set by default for CUDA PR testing")
set (EpetraExt_ENABLE_HDF5 OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Panzer_FADTYPE "Sacado::Fad::DFad<RealType>" CACHE STRING "Set by default for CUDA PR testing")
set (Kokkos_ENABLE_Debug_Bounds_Check ON CACHE BOOL "Set by default for CUDA PR testing")
set (KOKKOS_ENABLE_DEBUG ON CACHE BOOL "Set by default for CUDA PR testing")

# TPL settings specific to CUDA build
set (TPL_ENABLE_CUDA ON CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_CUSPARSE ON CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_Pthread OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_BinUtils OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_BLAS ON CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_BLAS_LIBRARIES "-L$ENV{BLAS_ROOT}/lib;-lopenblas;-lgfortran;-lgomp;-lm" CACHE STRING "Set by default for CUDA PR testing")
set (TPL_ENABLE_Boost ON CACHE BOOL "Set by default for CUDA PR testing")
set (Boost_INCLUDE_DIRS "$ENV{BOOST_ROOT}/include" CACHE FILEPATH "Set by default for CUDA PR testing")
set (TPL_Boost_LIBRARIES "$ENV{BOOST_ROOT}/lib/libboost_program_options.a;$ENV{BOOST_ROOT}/lib/libboost_system.a" CACHE FILEPATH  "Set by default for CUDA PR testing")
set (TPL_ENABLE_BoostLib ON CACHE BOOL "Set by default for CUDA PR testing")
set (BoostLib_INCLUDE_DIRS "$ENV{BOOST_ROOT}/include" CACHE FILEPATH "Set by default for CUDA PR testing")
set (TPL_BoostLib_LIBRARIES "$ENV{BOOST_ROOT}/lib/libboost_program_options.a;$ENV{BOOST_ROOT}/lib/libboost_system.a" CACHE FILEPATH "Set by default for CUDA PR testing")
set (TPL_ENABLE_METIS OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_ParMETIS OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_Scotch OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_HWLOC OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_LAPACK ON CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_LAPACK_LIBRARIES "-L$ENV{LAPACK_ROOT}/lib;-lopenblas;-lgfortran;-lgomp" CACHE STRING "Set by default for CUDA PR testing")
set (TPL_ENABLE_Zlib OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_CGNS OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_HDF5 ON CACHE BOOL "Set by default for CUDA PR testing")
set (HDF5_INCLUDE_DIRS "$ENV{HDF5_ROOT}/include" CACHE FILEPATH "Set by default for CUDA PR testing")
set (TPL_HDF5_LIBRARIES "-L$ENV{HDF5_ROOT}/lib;$ENV{HDF5_ROOT}/lib/libhdf5_hl.a;$ENV{HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl" CACHE FILEPATH "Set by default for CUDA PR testing")
set (TPL_ENABLE_Netcdf ON CACHE BOOL "Set by default for CUDA PR testing")
set (Netcdf_INCLUDE_DIRS "$ENV{NETCDF_C_INC}" CACHE FILEPATH "Set by default for CUDA PR testing")
IF ("$ENV{PARALLEL_NETCDF_ROOT}" STREQUAL "")
  SET(PARALLEL_NETCDF_ROOT "$ENV{NETCDF_C_ROOT}")
ELSE()
  SET(PARALLEL_NETCDF_ROOT "$ENV{PARALLEL_NETCDF_ROOT}")
ENDIF()
set (TPL_Netcdf_LIBRARIES "-L$ENV{NETCDF_C_ROOT}/lib;$ENV{NETCDF_C_ROOT}/lib/libnetcdf.a;${PARALLEL_NETCDF_ROOT}/lib/libpnetcdf.a;${TPL_HDF5_LIBRARIES}" CACHE STRING "Set by default for CUDA PR testing")
# SuperLU and SuperLUDist is available on ride and could be enabled for the CUDA PR build
set (TPL_ENABLE_SuperLU OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_SuperLUDist OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_BoostLib OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_Matio OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_DLlib_LIBRARIES "-ldl" CACHE FILEPATH "Set by default for CUDA PR testing")

# The compile times for two Panzer files went up to over 6 hours. This
# turns off one feature that allows these in about 24 minutes. Please remove
# when issue #7532 is resolved.
# Compile time issues addressed by #8377. Commenting out the override of 
# Sacado_NEW_FAD_DESIGN_IS_DEFAULT to return to default settings.
#set (Sacado_NEW_FAD_DESIGN_IS_DEFAULT OFF CACHE BOOL "Temporary fix for issue #7532" )

# Disable some packages that can't be tested with this PR build
set (Trilinos_ENABLE_ShyLU_NodeTacho OFF CACHE BOOL
  "Can't test Tacho with CUDA without RDC" FORCE)

# Force some tests to run in serial in this PR build (usually to resolve random timeouts that can occur under contention)
set (Intrepid2_unit-test_Discretization_Basis_HierarchicalBases_Hierarchical_Basis_Tests_MPI_1_SET_RUN_SERIAL ON CACHE BOOL "Run serial for CUDA PR testing")

# Options used to clean up build
set (MueLu_ParameterListInterpreterTpetra_MPI_1_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
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
#set (PanzerAdaptersSTK_CurlLaplacianExample-ConvTest-Quad-Order-4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
#set (PanzerAdaptersSTK_MixedPoissonExample-ConvTest-Hex-Order-3_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
# set (TrilinosCouplings_Example_Maxwell_MueLu_MPI_1_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
# set (TrilinosCouplings_Example_Maxwell_MueLu_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")

# Disable some tests that should not need to be disabled but do because the
# Trilinos PR tester is using too high a parallel level for ctest. (If the
# ctest parallel test level is dropped from 29 to 8, all of these will pass.)
set (MueLu_UnitTestsIntrepid2Tpetra_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
#set (PanzerAdaptersSTK_main_driver_energy-ss-blocked-tp_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
#set (PanzerAdaptersSTK_MixedCurlLaplacianExample-ConvTest-Tri-Order-1_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
#set (PanzerAdaptersSTK_PoissonInterfaceExample_3d_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
#set (PanzerMiniEM_MiniEM-BlockPrec_Augmentation_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
#set (PanzerMiniEM_MiniEM-BlockPrec_RefMaxwell_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ROL_example_PDE-OPT_poisson-boltzmann_example_01_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")

# Disable a couple of unit tests in test Kokkos_CoreUnitTest_Cuda_MPI_1 that
# are randomly failing in PR test iterations (#6799)
set (Kokkos_CoreUnitTest_Cuda_MPI_1_EXTRA_ARGS
  "--gtest_filter=-cuda.debug_pin_um_to_host:cuda.debug_serial_execution"
  CACHE STRING "Temporary disable for CUDA PR testing")


# Disable a few failing tests for initial release of Cuda 10.1.105 PR build
set (EpetraExt_inout_test_LL_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (EpetraExt_inout_test_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (Teko_testdriver_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (Zoltan2_fix4785_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")

set (Trilinos_ENABLE_COMPLEX ON CACHE BOOL "Testing for #9025")
set (Teuchos_ENABLE_COMPLEX ON CACHE BOOL "Testing for #9025")
set (Tpetra_INST_COMPLEX_DOUBLE ON CACHE BOOL "Testing for #9025")

# UVM = OFF
set (Kokkos_ENABLE_CUDA_UVM OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Tpetra_ENABLE_CUDA_UVM OFF CACHE BOOL "Set by default for CUDA PR testing")

# Turn off packages currently failing with UVM = OFF
set (Trilinos_ENABLE_Stokhos OFF CACHE BOOL "Turn off packages for non-UVM build")

# Turn off for sems-cuda-11 build
#set (Trilinos_ENABLE_Sacado OFF CACHE BOOL "Turn off packages for non-UVM build")
#set (Trilinos_ENABLE_SEACAS OFF CACHE BOOL "Turn off packages for non-UVM build")
#set (Trilinos_ENABLE_Ifpack2 OFF CACHE BOOL "Turn off packages for non-UVM build")
#set (Trilinos_ENABLE_Zoltan2 OFF CACHE BOOL "Turn off packages for non-UVM build")
#set (Trilinos_ENABLE_NOX OFF CACHE BOOL "Turn off packages for non-UVM build")
#set (Trilinos_ENABLE_Piro OFF CACHE BOOL "Turn off packages for non-UVM build")
#set (Trilinos_ENABLE_MueLu OFF CACHE BOOL "Turn off packages for non-UVM build")

# Turn off tests currently failing with UVM = OFF
# Packages with >5 failing tests
set (Anasazi_ENABLE_TESTS OFF CACHE BOOL "Turn off tests for non-UVM build")
set (Domi_ENABLE_TESTS OFF CACHE BOOL "Turn off tests for non-UVM build")
set (Kokkos_ENABLE_TESTS OFF CACHE BOOL "Turn off tests for non-UVM build")
set (KokkosKernels_ENABLE_TESTS OFF CACHE BOOL "Turn off tests for non-UVM build")
set (ROL_ENABLE_TESTS OFF CACHE BOOL "Turn off tests for non-UVM build")
set (SEACAS_ENABLE_TESTS OFF CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DD_ENABLE_TESTS OFF CACHE BOOL "Turn off tests for non-UVM build")
set (STK_ENABLE_TESTS OFF CACHE BOOL "Turn off tests for non-UVM build")


# ShyLU_DD UVM = OFF tests
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_one_rank_TLP_IPOU_DIM3_TPETRA_MPI_1_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_one_rank_TLP_GDSW_DIM2_TPETRA_MPI_1_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_one_rank_TLP_GDSW_DIM3_TPETRA_MPI_1_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_one_rank_TLP_RGDSW_DIM3_TPETRA_MPI_1_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_one_rank_TLP_RGDSW_DIM2_TPETRA_MPI_1_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_one_rank_TLP_IPOU_DIM2_TPETRA_MPI_1_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDBDDC_bddc_standard_interface_test_MPI_2_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDBDDC_bddc_simple_interface_test_MPI_2_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM2_DPN2_ORD0_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM2_DPN1_ORD1_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM2_DPN2_ORD1_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLBP_NB3_GDSW_DIM2_DPN1_ORD0_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM2_DPN1_ORD1_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM2_DPN2_ORD1_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_GDSW_DIM2_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM2_DPN1_ORD0_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_IPOUHarmonic_GDSW_DIM2_DPN2_ORD0_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_IPOUHarmonic_RGDSW_DIM2_DPN2_ORD1_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_IPOUHarmonic_RGDSW_DIM2_DPN2_ORD0_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLBP_NB2_GDSW_DIM2_DPN2_ORD0_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_IPOUHarmonic_GDSWStar_DIM2_DPN2_ORD0_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_GDSWStar_DIM2_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_GDSW_DIM2_DPN2_ORD0_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_GDSW_DIM2_DPN2_ORD1_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_elasticity_TLP_GDSW_DIM2_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_IPOUHarmonic_GDSWStar_DIM2_DPN2_ORD1_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_IPOUHarmonic_GDSWStar_DIM2_DPN1_ORD0_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_GDSWP_DIM2_DPN2_ORD1_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_RGDSW_DIM2_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_RGDSW_DIM2_DPN2_ORD1_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_RGDSW_DIM2_DPN2_ORD0_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_elasticity_GDSWP_DIM2_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_elasticity_TLP_RGDSW_DIM2_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_IPOUHarmonic_GDSW_DIM2_DPN2_ORD1_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_elasticity_RGDSWP_DIM2_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_RGDSW_DIM2_DPN1_ORD0_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_RGDSWP_DIM2_DPN2_ORD0_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_GDSWP_DIM2_DPN1_ORD0_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_GDSWP_DIM2_DPN2_ORD0_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_RGDSWP_DIM2_DPN1_ORD0_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_IPOUHarmonic_RGDSW_DIM2_DPN1_ORD0_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_IPOUHarmonic_GDSW_DIM2_DPN1_ORD0_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_RGDSWP_DIM2_DPN2_ORD1_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_GDSW_DIM2_DPN1_ORD0_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_overlap_Old_DIM2_DPN2_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_overlap_Matrix_DIM2_DPN2_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_overlap_Graph_DIM2_DPN2_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_overlap_Matrix_DIM2_DPN1_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_overlap_Graph_DIM2_DPN1_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_overlap_Old_DIM2_DPN1_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_interfacesets_DIM2_CSCommCrsGraph_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_interfacepartitionofunity_DIM2_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_interfacesets_DIM2_CSCommCrsMatrix_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_interfacesets_DIM2_CSOneToOne_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_GDSW_ML_DIM2_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_IPOUHarmonic_GDSW_ML_DIM2_DPN1_ORD0_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_RGDSW_ML_DIM2_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_elasticity_TLP_IPOUHarmonic_RGDSW_DIM2_TPETRA_DropCoupling_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDBDDC_bddc_sparse_solver_test_MPI_1_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_localpartitionofunitybasis_TPETRA_MPI_1_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_thyra_xpetra_laplace_TLP_GDSW_MUELU_DIM2_DPN2_ORD1_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_solverfactory_thyrasolver_belos_gmres_DIM2_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_solverfactory_thyrapreconditioner_frosch_onelevelpreconditioner_DIM2_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_solverfactory_muelu_DIM2_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_solverfactory_ifpack2_riluk_DIM2_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_solverfactory_froschpreconditioner_twolevelpreconditioner_DIM2_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_solverfactory_froschpreconditioner_twolevelblockpreconditioner_DIM2_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_solverfactory_belos_gmres_DIM2_TPETRA_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_solverfactory_amesos2_klu_DIM3_TPETRA_MPI_1_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ShyLU_DDFROSch_test_solverfactory_amesos2_klu_DIM2_TPETRA_MPI_1_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")

# Adelus UVM = OFF tests
set (Adelus_vector_random_MPI_1_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (Adelus_vector_random_MPI_2_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (Adelus_vector_random_MPI_3_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (Adelus_vector_random_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
# Thyra UVM = OFF tests
set (ThyraTpetraAdapters_Simple2DTpetraModelEvaluatorUnitTests_MPI_1_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ThyraTpetraAdapters_TpetraThyraWrappersUnitTests_serial_MPI_1_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (ThyraTpetraAdapters_TpetraThyraWrappersUnitTests_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
# Piro UVM = OFF tests
set (Piro_AnalysisDriverTpetra_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (Piro_ThyraSolverTpetra_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
# MiniTensor UVM = OFF tests
set (MiniTensor_test_Test_01_MPI_1_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (MiniTensor_test_Test_02_MPI_1_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
# SEACAS UVM = OFF tests
set (SEACASIoss_exodus32_to_exodus32_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (SEACASIoss_exodus_fpp_serialize_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
# Stratimikos UVM = OFF tests
set (Stratimikos_test_single_amesos2_tpetra_solver_driver_KLU2_MPI_1_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")

# MueLu UVM = OFF tests
set (MueLu_ReitzingerPFactory_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")
set (MueLu_Maxwell3D-Tpetra_2_MPI_4_DISABLE ON CACHE BOOL "Turn off tests for non-UVM build")

include("${CMAKE_CURRENT_LIST_DIR}/PullRequestLinuxCommonTestingSettings.cmake")
