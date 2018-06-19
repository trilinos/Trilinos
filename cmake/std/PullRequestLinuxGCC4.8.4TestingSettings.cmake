# This file contains the options needed to both run the pull request testing
# for Trilinos for the Linux GCC 4.8.4 pull request testing builds, and to reproduce
# the errors reported by those builds. Prior to using this this file, the
# appropriate set of SEMS modules must be loaded and accessible through the
# SEMS NFS mount. (See the sems/PullRequestGCC*TestingEnv.sh files.)

# Usage: cmake -C PullRequestLinuxGCC4.8.4TestingSettings.cmake 

# Misc options typically added by CI testing mode in TriBITS

# Use the below option only when submitting to the dashboard
#set (CTEST_USE_LAUNCHERS ON CACHE BOOL "Set by default for PR testing")

set (Trilinos_ENABLE_OpenMP ON CACHE BOOL "Set by default for PR testing")
set (MPI_EXEC_PRE_NUMPROCS_FLAGS "--bind-to;none" CACHE STRING "Set by default for PR testing")
# NOTE: The above is a workaround for the problem of having threads on MPI
# ranks bind to the same cores (see #2422).
 
# Disable just one Teko sub-unit test that fails with GCC 4.8.4 + OpenMP (#2712)
set (Teko_DISABLE_LSCSTABALIZED_TPETRA_ALPAH_INV_D ON CACHE BOOL "Temporarily disabled in PR testing")

include("${CMAKE_CURRENT_LIST_DIR}/PullRequestLinuxGCCTestingSettings.cmake")

