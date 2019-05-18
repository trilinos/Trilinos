# This file contains the options needed to both run the pull request testing
# for Trilinos for the Linux GCC 7.2.0 pull request testing builds, and to reproduce
# the errors reported by those builds. Prior to using this this file, the
# appropriate set of SEMS modules must be loaded and accessible through the
# SEMS NFS mount. (See the sems/PullRequestGCC*TestingEnv.sh files.)

# Usage: cmake -C PullRequestLinuxGCC7.2.0TestingSettings.cmake

# Misc options typically added by CI testing mode in TriBITS

# Use the below option only when submitting to the dashboard
#set (CTEST_USE_LAUNCHERS ON CACHE BOOL "Set by default for PR testing")

set (MPI_EXEC_PRE_NUMPROCS_FLAGS "--bind-to;none" CACHE STRING "Set by default for PR testing")
# NOTE: The above is a workaround for the problem of having threads on MPI
# ranks bind to the same cores (see #2422).

set (Trilinos_ENABLE_COMPLEX_DOUBLE ON CACHE BOOL "Set by default for PR testing to exercise complex doubles case")

# Disable just one Teko sub-unit test that fails with openmpi 1.10 (#2712)
set (Teko_DISABLE_LSCSTABALIZED_TPETRA_ALPAH_INV_D ON CACHE BOOL "Temporarily disabled in PR testing")

include("${CMAKE_CURRENT_LIST_DIR}/PullRequestLinuxCommonTestingSettings.cmake")

# Adding warnings as errors flags to this PR build
set(CMAKE_CXX_FLAGS "-Wall -Wno-clobbered -Wno-vla -Wno-pragmas -Wno-unknown-pragmas -Wno-unused-local-typedefs -Wno-literal-suffix -Wno-deprecated-declarations -Wno-misleading-indentation -Wno-int-in-bool-context -Wno-maybe-uninitialized -Wno-nonnull-compare -Wno-address -Wno-inline -DTRILINOS_HIDE_DEPRECATED_HEADER_WARNINGS" CACHE STRING "Warning settings")

#set (Amesos_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Amesos2_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Anasazi_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (AztecOO_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Belos_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Domi_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (Epetra_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (EpetraExt_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (FEI_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Galeri_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (GlobiPack_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Gtest_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Ifpack_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Ifpack2_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Intrepid_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Intrepid2_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Isorropia_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Kokkos_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (KokkosKernels_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (MiniTensor_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (ML_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (MueLu_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (NOX_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (OptiPack_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Pamgen_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Panzer_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Phalanx_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Pike_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Piro_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
#set (ROL_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (RTOp_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Rythmos_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Sacado_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (SEACAS_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Shards_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (ShyLU_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (ShyLU_DD_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (ShyLU_Node_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (STK_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Stokhos_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Stratimikos_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Teko_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Tempus_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Teuchos_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (ThreadPool_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Thyra_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Tpetra_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (TrilinosCouplings_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (TrilinosSS_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Triutils_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Xpetra_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Zoltan_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
set (Zoltan2_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" CACHE STRING "Warnings as errors setting")
