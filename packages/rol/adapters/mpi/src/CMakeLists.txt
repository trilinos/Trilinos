

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

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# All headers
#

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/vector)
  APPEND_SET(HEADERS 
    vector/ROL_PinTVector.hpp
    vector/ROL_PinTCommunicators.hpp
    vector/ROL_PinTVectorCommunication.hpp
    vector/ROL_PinTVectorCommunication_StdVector.hpp
)

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/function)
  APPEND_SET(HEADERS 
    function/ROL_PinTConstraint.hpp
)









INSTALL(FILES ${HEADERS} DESTINATION ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})
