# This file contains the options needed to both run the pull request testing
# for Trilinos for the Linux Intel 19 pull request testing builds, and to reproduce
# the errors reported by those builds. Prior to using this this file, the
# appropriate set of SEMS modules must be loaded and accessible through the
# SEMS NFS mount. (See the sems/PullRequestGCC*TestingEnv.sh files.)

# Usage: cmake -C PullRequestLinuxIntel19.0.5TestingSettings.cmake

set (CMAKE_CXX_STANDARD "14" CACHE STRING "Set C++ standard to C++14")
#Failing tests under C++14
#Remove line if test has been fixed
set (Piro_AnalysisDriver_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (Piro_AnalysisDriverTpetra_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_adapters_epetra_test_sol_EpetraSROMSampleGenerator_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_adapters_minitensor_test_function_test_01_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_adapters_minitensor_test_sol_test_01_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_adapters_teuchos_test_sol_solSROMGenerator_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_example_burgers-control_example_01_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_example_diode-circuit_example_01_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_example_parabolic-control_example_01_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_example_parabolic-control_example_02_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_example_parabolic-control_example_03_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_example_parabolic-control_example_04_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_example_parabolic-control_example_05_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_example_PDE-OPT_0ld_adv-diff-react_example_02_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_example_PDE-OPT_0ld_poisson_example_01_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_example_PDE-OPT_0ld_stoch-adv-diff_example_01_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_example_PDE-OPT_adv-diff-react_example_02_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_example_PDE-OPT_navier-stokes_example_02_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_example_PDE-OPT_nonlinear-elliptic_example_01_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_example_PDE-OPT_nonlinear-elliptic_example_02_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_example_PDE-OPT_obstacle_example_01_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_example_PDE-OPT_topo-opt_poisson_example_01_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_example_poisson-control_example_01_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_example_poisson-control_example_02_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_example_poisson-inversion_example_02_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_example_tensor-opt_example_01_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_NonlinearProblemTest_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_test_algorithm_OptimizationSolverStatusTestInput_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_test_elementwise_BoundConstraint_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_test_function_BinaryConstraintCheck_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_test_function_ExplicitLinearConstraintCheck_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_test_step_AugmentedLagrangianStep_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_test_step_BoxConstrained_LineSearch_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_test_step_BoxConstrained_LM_TrustRegion_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_test_step_BoxConstrained_PrimalDualActiveSet_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_test_step_BoxConstrained_TrustRegion_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_test_step_CubicTest_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_test_step_fletcher_ALLPROBLEMS_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_test_step_fletcher_BOUNDFLETCHER_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_test_step_FletcherStep_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_test_step_interiorpoint_PrimalDualNewtonKrylov_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_test_step_InteriorPointStep_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_test_step_LineSearch_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_test_step_MoreauYosidaPenaltyStep_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_test_step_TrustRegion_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (ROL_tutorial_BoundAndInequality_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")


# Tests turned off until an issue with Scotch can be resolved
set (Zoltan2_scotch_example_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (Zoltan2_Partitioning1_VWeights_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (Zoltan2_Partitioning1_OneProc_VWeights_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (Zoltan2_Partitioning1_OneProc_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (Zoltan2_Partitioning1_OneProc_EWeights_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (Zoltan2_Partitioning1_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (Zoltan2_Partitioning1_EWeights_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (Zoltan2_pamgenMeshAdapterTest_scotch_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (Phalanx_dynamic_data_layout_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (Zoltan_ENABLE_Scotch OFF CACHE BOOL "Temporarily disabled in PR testing")

# Misc options typically added by CI testing mode in TriBITS

# Use the below option only when submitting to the dashboard
#set (CTEST_USE_LAUNCHERS ON CACHE BOOL "Set by default for PR testing")

set (MueLu_UnitTestsEpetra_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (MueLu_UnitTestsEpetra_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (MueLu_UnitTestsTpetra_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")

# (Temporarily) Disable randomly failing ROL test (#3103)
set (ROL_example_poisson-inversion_example_01_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")

include("${CMAKE_CURRENT_LIST_DIR}/PullRequestLinuxCommonTestingSettings.cmake")

set (Tpetra_INST_INT_INT ON CACHE BOOL "INST_INT_INT ON")
set (Trilinos_ENABLE_STKBalance OFF CACHE BOOL "Hard disabled since Tpetra_INST_INT_INT=ON in this build" FORCE)
#STK-TODO: try to remember to come back and remove this when stk-balance
#is able to tolerate int as a global-index.

set (TPL_Netcdf_LIBRARIES "-L${Netcdf_LIBRARY_DIRS}/lib;${Netcdf_LIBRARY_DIRS}/libnetcdf.so;${Netcdf_LIBRARY_DIRS}/libpnetcdf.a" CACHE STRING "Set by default for CUDA PR testing")

set(CMAKE_CXX_FLAGS "-fPIC -Wall -Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat -Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses -Wpointer-sign -Wreorder -Wreturn-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable -Wwrite-strings" CACHE STRING "enable relocatable code, turn on WError settings")

# -lifcore is needed for CMake > 3.10.x builds.
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lifcore" CACHE STRING "updated by Pull Request")

#set (Anasazi_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Belos_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Domi_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Epetra_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (EpetraExt_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (FEI_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Ifpack_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Ifpack2_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Intrepid2_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Kokkos_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (KokkosKernels_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (ML_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (MueLu_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (NOX_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Panzer_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Phalanx_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Pike_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Piro_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (ROL_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Sacado_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (SEACAS_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Shards_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Stokhos_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Tempus_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Teuchos_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Tpetra_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Triutils_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Zoltan2_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
