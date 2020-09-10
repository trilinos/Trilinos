# This file contains the options needed to both run the pull request testing
# for Trilinos for the Linux GCC 4.8.4 pull request testing builds, and to reproduce
# the errors reported by those builds. Prior to using this this file, the
# appropriate set of SEMS modules must be loaded and accessible through the
# SEMS NFS mount. (See the sems/PullRequestGCC*TestingEnv.sh files.)

# Usage: cmake -C PullRequestLinuxGCC4.8.4TestingSettings.cmake

# Misc options typically added by CI testing mode in TriBITS

# Use the below option only when submitting to the dashboard
#set (CTEST_USE_LAUNCHERS ON CACHE BOOL "Set by default for PR testing")

set (TFW_Python3_Testing ON CACHE BOOL "Set by default for PR testing")
set (TFW_Python_Testing ON CACHE BOOL "Set by default for PR testing")

#set(PYTHON_EXECUTABLE /projects/sierra/linux_rh7/install/Python/3.6.3/bin/python CACHE FILEPATH "Set by default for PR testing")
#set(PYTHON_EXECUTABLE /projects/sierra/linux_rh7/install/Python/3.6.10/bin/python CACHE FILEPATH "Set by default for PR testing")
set(PYTHON_PIP_EXECUTABLE "pip3" CACHE STRING "Set by default for PR testing")

set(PYTHON_EXECUTABLE_SEARCH_PATHS
    /projects/sierra/linux_rh7/install/Python/3.6.3
    /projects/sierra/linux_rh7/install/Python/3.6.10
)
find_program(PYTHON_EXECUTABLE
    NAMES python3 python
    PATHS ${PYTHON_EXECUTABLE_SEARCH_PATHS}
    PATH_SUFFIXES bin
    DOC "Set by default for PR testing"
    NO_DEFAULT_PATH
    )
if(DEFINED PYTHON_EXECUTABLE_NOTFOUND)
    message(FATAL_ERROR "Unable to locate Python in ${PYTHON_EXECUTABLE_SEARCH_PATHS}")
else()
    message(STATUS "PYTHON FOUND: ${PYTHON_EXECUTABLE}")
endif()


set (TPL_ENABLE_Boost OFF CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_BoostLib OFF CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_ParMETIS OFF CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_Zlib OFF CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_HDF5 OFF CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_Netcdf OFF CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_SuperLU OFF CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_Scotch OFF CACHE BOOL "Set by default for PR testing")

include("${CMAKE_CURRENT_LIST_DIR}/PullRequestLinuxCommonTestingSettings.cmake")

