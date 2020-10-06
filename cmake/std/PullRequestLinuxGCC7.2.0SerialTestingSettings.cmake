# This file contains the options needed to both run the pull request testing
# for Trilinos for the Linux GCC 7.2.0 pull request testing builds, and to reproduce
# the errors reported by those builds. Prior to using this this file, the
# appropriate set of SEMS modules must be loaded and accessible through the
# SEMS NFS mount. (See the sems/PullRequestGCC*TestingEnv.sh files.)

# Usage: cmake -C PullRequestLinuxGCC7.2.0TestingSettings.cmake

set (CMAKE_CXX_STANDARD "14" CACHE STRING "Set C++ standard to C++14")

# Misc options typically added by CI testing mode in TriBITS

# Use the below option only when submitting to the dashboard
#set (CTEST_USE_LAUNCHERS ON CACHE BOOL "Set by default for PR testing")

# Disabling fortran for this build as that is a customer use case
set (Trilinos_ENABLE_Fortran OFF CACHE BOOL "Set by default for SERIAL PR test")

set (Trilinos_ENABLE_COMPLEX_DOUBLE ON CACHE BOOL "Set by default for PR testing to exercise complex doubles case")

include("${CMAKE_CURRENT_LIST_DIR}/PullRequestLinuxCommonTestingSettingsSERIAL.cmake")

# Adding warnings as errors flags to this PR build
set(CMAKE_CXX_FLAGS "-Wall -Wno-clobbered -Wno-vla -Wno-pragmas -Wno-unknown-pragmas -Wno-unused-local-typedefs -Wno-literal-suffix -Wno-deprecated-declarations -Wno-misleading-indentation -Wno-int-in-bool-context -Wno-maybe-uninitialized -Wno-nonnull-compare -Wno-address -Wno-inline -Wno-unused-but-set-variable -Wno-unused-variable -Wno-unused-label -Werror -DTRILINOS_HIDE_DEPRECATED_HEADER_WARNINGS" CACHE STRING "Warnings as errors settings")
