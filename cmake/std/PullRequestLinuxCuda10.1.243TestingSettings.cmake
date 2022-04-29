# This file contains the options needed to both run the pull request testing
# for Trilinos for the CUDA 10.1.105 pull request testing builds, and to reproduce
# the errors reported by those builds. Prior to using this this file, the
# appropriate set of modules must be loaded and path must be augmented.
# (See the sems/PullRequestCuda10.1.105TestingEnv.sh files.)

# Usage: cmake -C PullRequestLinuxCUDA10.1.243TestingSettings.cmake

# Misc options typically added by CI testing mode in TriBITS

# Use the below option only when submitting to the dashboard
set (CTEST_USE_LAUNCHERS ON CACHE BOOL "Set by default for PR testing")
set (MPI_EXEC "${CMAKE_CURRENT_LIST_DIR}/sems/trilinos_jsrun" CACHE STRING "We use this wrapper to get the cuda hooks turned off") 

# Options necessary for CUDA build
set (TPL_ENABLE_MPI ON CACHE BOOL "Set by default for CUDA PR testing")
set (Kokkos_ENABLE_CUDA ON CACHE BOOL "Set by default for CUDA PR testing")
set (Kokkos_ENABLE_CUDA_UVM ON CACHE BOOL "Set by default for CUDA PR testing")
set (Kokkos_ARCH_POWER9 ON CACHE BOOL "Set by default for CUDA PR testing")
set (Kokkos_ARCH_VOLTA70 ON CACHE BOOL "Set by default for CUDA PR testing")
set (Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Sacado_ENABLE_HIERARCHICAL_DFAD ON CACHE BOOL "Set by default for CUDA PR testing")
set (Kokkos_ENABLE_CXX11_DISPATCH_LAMBDA ON CACHE BOOL "Set by default for CUDA PR testing")
set (Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "Set by default for CUDA PR testing")
set (Phalanx_KOKKOS_DEVICE_TYPE CUDA CACHE STRING "Set by default for CUDA PR testing")
set (MPI_EXEC_POST_NUMPROCS_FLAGS "--rs_per_socket;4" CACHE STRING "Set by default for CUDA PR testing")
set (MPI_EXEC_NUMPROCS_FLAG "-p" CACHE STRING "Set by default for CUDA PR testing")

# Options set to match the ATDM build
set (Trilinos_ENABLE_DEBUG OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Trilinos_ENABLE_DEBUG_SYMBOLS OFF CACHE BOOL "Set by default for CUDA PR testing")
set (BUILD_SHARED_LIBS OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Tpetra_INST_SERIAL ON CACHE BOOL "Set by default for CUDA PR testing")
## this was moved to the environment settings in sems/PullRequestCuda10.1.243TestingEnv.sh
## set (TPETRA_ASSUME_CUDA_AWARE_MPI ON CACHE BOOL "Set by default for CUDA PR testing")

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
set (TPL_BLAS_LIBRARY_DIRS "$ENV{CBLAS_ROOT}" "Set by default for CUDA PR testing")
set (TPL_BLAS_LIBRARIES "-L$ENV{CBLAS_ROOT}/lib;-lblas;-lgfortran;-lgomp;-lm" CACHE STRING "Set by default for CUDA PR testing")
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
set (TPL_LAPACK_LIBRARIES "-L$ENV{LAPACK_ROOT}/lib;-llapack;-lgfortran;-lgomp" CACHE STRING "Set by default for CUDA PR testing")
set (TPL_ENABLE_Zlib OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_CGNS OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_HDF5 ON CACHE BOOL "Set by default for CUDA PR testing")
set (HDF5_INCLUDE_DIRS "$ENV{HDF5_ROOT}/include" CACHE FILEPATH "Set by default for CUDA PR testing")
set (TPL_HDF5_LIBRARIES "-L$ENV{HDF5_ROOT}/lib;$ENV{HDF5_ROOT}/lib/libhdf5_hl.a;$ENV{HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl" CACHE FILEPATH "Set by default for CUDA PR testing")
set (TPL_ENABLE_Netcdf ON CACHE BOOL "Set by default for CUDA PR testing")
set (Netcdf_INCLUDE_DIRS "$ENV{NETCDF_ROOT}/include" CACHE FILEPATH "Set by default for CUDA PR testing")
IF ("$ENV{PNETCDF_ROOT}" STREQUAL "")
  SET(PNETCDF_ROOT "$ENV{NETCDF_ROOT}")
ELSE()
  SET(PNETCDF_ROOT "$ENV{PNETCDF_ROOT}")
ENDIF()
set (TPL_Netcdf_LIBRARIES "-L$ENV{NETCDF_ROOT}/lib64;$ENV{NETCDF_ROOT}/lib64/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${TPL_HDF5_LIBRARIES}" CACHE STRING "Set by default for CUDA PR testing")

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

# Temporary options to clean up build
set (Teko_ModALPreconditioner_MPI_1_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (MueLu_ParameterListInterpreterTpetra_MPI_1_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (MueLu_ParameterListInterpreterTpetraHeavy_MPI_1_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (STKUnit_tests_stk_ngp_test_utest_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
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
set (ROL_NonlinearProblemTest_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (PanzerAdaptersSTK_CurlLaplacianExample-ConvTest-Quad-Order-4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (PanzerAdaptersSTK_MixedPoissonExample-ConvTest-Hex-Order-3_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (TrilinosCouplings_Example_Maxwell_MueLu_MPI_1_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (TrilinosCouplings_Example_Maxwell_MueLu_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")

# Disable some tests that should not need to be disabled but do because the
# Trilinos PR tester is using too high a parallel level for ctest. (If the
# ctest parallel test level is dropped from 29 to 8, all of these will pass.)
set (MueLu_UnitTestsIntrepid2Tpetra_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (PanzerAdaptersSTK_main_driver_energy-ss-blocked-tp_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (PanzerAdaptersSTK_MixedCurlLaplacianExample-ConvTest-Tri-Order-1_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (PanzerAdaptersSTK_PoissonInterfaceExample_3d_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (PanzerMiniEM_MiniEM-BlockPrec_Augmentation_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (PanzerMiniEM_MiniEM-BlockPrec_RefMaxwell_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ROL_example_PDE-OPT_poisson-boltzmann_example_01_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")

# Disable a couple of unit tests in test KokkosCore_UnitTest_Cuda_MPI_1 that
# are randomly failing in PR test iterations (#6799)
set (KokkosCore_UnitTest_Cuda_MPI_1_EXTRA_ARGS
  "--gtest_filter=-cuda.debug_pin_um_to_host:cuda.debug_serial_execution"
  CACHE STRING "Temporary disable for CUDA PR testing")

# Disable a few failing tests for initial release of Cuda 10.1.243 PR build
# these fail due to a bug in spectrum-mpi : see issue #8798 for details
set (Zoltan_ch_simple_zoltan_parallel_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (Zoltan_hg_simple_zoltan_parallel_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (Zoltan_mpiMinLoc_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")

# Disable SEACAS tests that grep results out of stderr...
set (SEACASIoss_create_path_fpp_DISABLE ON CACHE BOOL "Temporary disable due to jsrun polluting stderr")
set (SEACASAprepro_lib_aprepro_lib_unit_test_DISABLE ON CACHE BOOL "Temporary disable due to jsrun polluting stderr")
set (SEACASAprepro_lib_aprepro_lib_array_test_DISABLE ON CACHE BOOL "Temporary disable due to jsrun polluting stderr")
set (SEACASAprepro_aprepro_unit_test_DISABLE ON CACHE BOOL "Temporary disable due to jsrun polluting stderr")
set (SEACASAprepro_aprepro_array_test_DISABLE ON CACHE BOOL "Temporary disable due to jsrun polluting stderr")
set (SEACASAprepro_aprepro_command_line_vars_test_DISABLE ON CACHE BOOL "Temporary disable due to jsrun polluting stderr")
set (SEACASAprepro_aprepro_command_line_include_test_DISABLE ON CACHE BOOL "Temporary disable due to jsrun polluting stderr")
set (SEACASAprepro_aprepro_test_dump_reread_DISABLE ON CACHE BOOL "Temporary disable due to jsrun polluting stderr")

# Disable tests detailed in issue #8796
set (Adelus_vector_random_MPI_1_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")

# Disable tests detailed in issue #8799
set (ROL_adapters_minitensor_test_function_test_01_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ROL_adapters_minitensor_test_function_test_02_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ROL_adapters_minitensor_test_sol_test_01_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ROL_adapters_minitensor_test_vector_test_01_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (ROL_test_algorithm_TypeB_TrustRegionSPG_MPI_1_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")

# Disable tests detailed in issue #8800
set (KokkosCore_UnitTest_CudaTimingBased_MPI_1_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")

# Disable tests detailed in issue #8129
set (PanzerDiscFE_integration_values2_MPI_1_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")
set (Stokhos_TpetraCrsMatrixMPVectorUnitTest_Cuda_MPI_4_DISABLE ON CACHE BOOL "Temporary disable for CUDA PR testing")

set (CMAKE_CXX_STANDARD "14" CACHE STRING "Set C++ standard to C++14")
# set (CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Kokkos turns off CXX extensions")

include("${CMAKE_CURRENT_LIST_DIR}/PullRequestLinuxCommonTestingSettings.cmake")
