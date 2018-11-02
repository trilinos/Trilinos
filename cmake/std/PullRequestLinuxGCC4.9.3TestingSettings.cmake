# This file contains the options needed to both run the pull request testing
# for Trilinos for the Linux GCC 4.9.3 pull request testing builds, and to reproduce
# the errors reported by those builds. Prior to using this this file, the
# appropriate set of SEMS modules must be loaded and accessible through the
# SEMS NFS mount. (See the sems/PullRequestGCC*TestingEnv.sh files.)

# Usage: cmake -C PullRequestLinuxGCC4.9.3TestingSettings.cmake

# Misc options typically added by CI testing mode in TriBITS

# Use the below option only when submitting to the dashboard
#set (CTEST_USE_LAUNCHERS ON CACHE BOOL "Set by default for PR testing")

set (MPI_EXEC_PRE_NUMPROCS_FLAGS "--bind-to;none" CACHE STRING "Set by default for PR testing")
# NOTE: The above is a workaround for the problem of having threads on MPI
# ranks bind to the same cores (see #2422).

# Disable just one Teko sub-unit test that fails with openmpi 1.10 (#2712)
set (Teko_DISABLE_LSCSTABALIZED_TPETRA_ALPAH_INV_D ON CACHE BOOL "Temporarily disabled in PR testing")

include("${CMAKE_CURRENT_LIST_DIR}/PullRequestLinuxCommonTestingSettings.cmake")

#Adding warnings as errors flags to this PR build
#This should fail. Starting to break down package by package.
set (CMAKE_CXX_FLAGS "-Wall -ansi -pedantic -Werror -Wno-unknown-pragmas -Wno-narrowing -Wno-pragmas -Wno-delete-non-virtual-dtor" CACHE STRING "Warnings as errors setting")
#set (Teuchos_CXX_FLAGS "-Wall -ansi -pedantic -Werror -Wno-unknown-pragmas -Wno-narrowing -Wno-pragmas -Wno-delete-non-virtual-dtor" CACHE STRING "Warnings as errors setting")
#set (Intrepid2_CXX_FLAGS "-Wall -ansi -pedantic -Werror -Wno-unknown-pragmas -Wno-narrowing -Wno-pragmas -Wno-delete-non-virtual-dtor" CACHE STRING "Warnings as errors setting")
#set (KokkosKernels_CXX_FLAGS "-Wall -ansi -pedantic -Werror -Wno-unknown-pragmas -Wno-narrowing -Wno-pragmas -Wno-delete-non-virtual-dtor" CACHE STRING "Warnings as errors setting")
#set (Phalanx_CXX_FLAGS "-Wall -ansi -pedantic -Werror -Wno-unknown-pragmas -Wno-narrowing -Wno-pragmas -Wno-delete-non-virtual-dtor" CACHE STRING "Warnings as errors setting")
#set (Pike_CXX_FLAGS "-Wall -ansi -pedantic -Werror -Wno-unknown-pragmas -Wno-narrowing -Wno-pragmas -Wno-delete-non-virtual-dtor" CACHE STRING "Warnings as errors setting")
#set (ShyLU_Node_CXX_FLAGS "-Wall -ansi -pedantic -Werror -Wno-unknown-pragmas -Wno-narrowing -Wno-pragmas -Wno-delete-non-virtual-dtor" CACHE STRING "Warnings as errors setting")
#set (Shards_CXX_FLAGS "-ansi -Werror -Wno-unknown-pragmas -Wno-narrowing -Wno-pragmas -Wno-delete-non-virtual-dtor" CACHE STRING "Warnings as errors setting")
#set (Zoltan_CXX_FLAGS "-Wall -ansi -pedantic -Werror -Wno-unknown-pragmas -Wno-narrowing -Wno-pragmas -Wno-delete-non-virtual-dtor" CACHE STRING "Warnings as errors setting")
#set (Epetra_CXX_FLAGS "-ansi -Werror -Wno-unknown-pragmas -Wno-narrowing -Wno-pragmas -Wno-delete-non-virtual-dtor" CACHE STRING "Warnings as errors setting")
