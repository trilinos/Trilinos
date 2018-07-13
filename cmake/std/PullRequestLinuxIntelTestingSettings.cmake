# This file contains the options needed to both run the pull request testing
# for Trilinos for the Linux Intel 17 pull request testing builds, and to reproduce
# the errors reported by those builds. Prior to using this this file, the
# appropriate set of SEMS modules must be loaded and accessible through the
# SEMS NFS mount. (See the sems/PullRequestGCC*TestingEnv.sh files.)

# Usage: cmake -C PullRequestLinuxIntelTestingSettings.cmake

# Misc options typically added by CI testing mode in TriBITS

# Use the below option only when submitting to the dashboard
#set (CTEST_USE_LAUNCHERS ON CACHE BOOL "Set by default for PR testing")

set (MueLu_UnitTestsEpetra_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (MueLu_UnitTestsEpetra_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (MueLu_UnitTestsTpetra_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (MueLu_UnitTestsTpetra_MPI_4_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (KokkosCore_UnitTest_Serial_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")
set (Zoltan2_simplePamgenTest_MPI_3_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")

# (Temporarily) Disable randomly failing ROL test (#3103)
set (ROL_example_poisson-inversion_example_01_MPI_1_DISABLE ON CACHE BOOL "Temporarily disabled in PR testing")

include("${CMAKE_CURRENT_LIST_DIR}/PullRequestLinuxCommonTestingSettings.cmake")
