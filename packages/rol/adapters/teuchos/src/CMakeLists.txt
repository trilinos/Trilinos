

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
    vector/ROL_CArrayVector.hpp
    vector/ROL_TeuchosVector.hpp
    vector/ROL_MultiVector.hpp
    vector/ROL_MultiVectorDefault.hpp
)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/function)
  APPEND_SET(HEADERS 
    function/ROL_Teuchos_Objective.hpp
    function/ROL_Teuchos_Constraint.hpp
)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sol)
  APPEND_SET(HEADERS 
    sol/ROL_SingletonTeuchosBatchManager.hpp
    sol/ROL_StdTeuchosBatchManager.hpp
    sol/ROL_TeuchosBatchManager.hpp
)




INSTALL(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})
