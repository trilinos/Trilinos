###############################################################################
#
#              ATDM Configuration for all Trilinos packages
#
###############################################################################

INCLUDE("${CMAKE_CURRENT_LIST_DIR}/utils/ATDMDevEnvUtils.cmake")

# Only enable Primary Tested packages for this build
ATDM_SET_ENABLE(${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE OFF)

# Need Fortran to build some SEACAS code
ATDM_SET_ENABLE(${PROJECT_NAME}_ENABLE_Fortran ON)

# Don't have stuff for Matio and some SEACAS subpackage has required
# dependency on this
ATDM_SET_ENABLE(TPL_ENABLE_Matio OFF)

# Now include the rest of the settings
INCLUDE("${CMAKE_CURRENT_LIST_DIR}/ATDMDevEnvSettings.cmake")
