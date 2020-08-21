# This file contains the options needed to both run the pull request testing
# for Trilinos for the Linux GCC 8.3.0 pull request testing builds, and to reproduce
# the errors reported by those builds. Prior to using this this file, the
# appropriate set of SEMS modules must be loaded and accessible through the
# SEMS NFS mount. (See the sems/PullRequestGCC*TestingEnv.sh files.)

# Usage: cmake -C PullRequestLinuxGCC8.3.0TestingSettings.cmake

# Misc options typically added by CI testing mode in TriBITS

# Use the below option only when submitting to the dashboard
#set (CTEST_USE_LAUNCHERS ON CACHE BOOL "Set by default for PR testing")

set (MPI_EXEC_PRE_NUMPROCS_FLAGS "--bind-to;none" CACHE STRING "Set by default for PR testing")
# NOTE: The above is a workaround for the problem of having threads on MPI
# ranks bind to the same cores (see #2422).

set (Trilinos_ENABLE_COMPLEX_DOUBLE ON CACHE BOOL "Set by default for PR testing to exercise complex doubles case")

# Enable C++14 for this build. In the future how we do this will change, but for
# now we need to use the Kokkos option

set (CMAKE_CXX_STANDARD "14" CACHE STRING "Set C++ standard to C++14")
set (CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Kokkos turns off CXX extensions")

# Disable just one Teko sub-unit test that fails with openmpi 1.10 (#2712)
set (Teko_DISABLE_LSCSTABALIZED_TPETRA_ALPAH_INV_D ON CACHE BOOL "Temporarily disabled in PR testing")

include("${CMAKE_CURRENT_LIST_DIR}/PullRequestLinuxCommonTestingSettings.cmake")

set(CMAKE_CXX_FLAGS "-std=c++14 -fno-strict-aliasing -Wall -Wno-clobbered -Wno-vla -Wno-pragmas -Wno-unknown-pragmas -Wno-parentheses -Wno-unused-local-typedefs -Wno-literal-suffix -Wno-deprecated-declarations -Wno-misleading-indentation -Wno-int-in-bool-context -Wno-maybe-uninitialized -Wno-class-memaccess -Wno-inline -Wno-nonnull-compare -Wno-address -DTRILINOS_HIDE_DEPRECATED_HEADER_WARNINGS" CACHE STRING "Warnings as errors settings")



