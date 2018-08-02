

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

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/vector)
  APPEND_SET(HEADERS 
    vector/ROL_PinTVector.hpp
)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/function)
  APPEND_SET(HEADERS 
    function/ODEConstraint_TimeSimOpt.hpp
    function/ROL_Constraint_PinTSimOpt.hpp
    function/ROL_PinTConstraint.hpp
)









INSTALL(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})
