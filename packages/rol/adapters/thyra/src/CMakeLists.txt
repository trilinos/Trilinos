
INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

# TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# All headers
#

#
# vector headers
#

IF (ROL_ENABLE_Thyra)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/vector)
  APPEND_SET(HEADERS
    vector/ROL_ThyraVector.hpp
  )
ENDIF()

#
# function headers
#

IF (ROL_ENABLE_Thyra)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/function)
  APPEND_SET(HEADERS
    function/ROL_ThyraME_Objective.hpp
    function/ROL_Thyra_BoundConstraint.hpp
    function/TOpEleWiseBoundsHelpers.hpp
  )
ENDIF()

INSTALL(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})
