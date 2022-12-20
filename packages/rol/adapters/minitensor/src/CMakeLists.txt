

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

#
# vector headers
#

IF (ROL_ENABLE_Boost AND ROL_ENABLE_MiniTensor)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/vector)
  APPEND_SET(HEADERS
    vector/ROL_MiniTensor_Vector.hpp
    vector/ROL_MiniTensor_Vector_Def.hpp
  )
ENDIF()

#
# function headers
#

IF (ROL_ENABLE_Boost AND ROL_ENABLE_MiniTensor)
#  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/function)
#  APPEND_SET(HEADERS
#    function/ROL_MiniTensor_BoundConstraint.hpp
#    function/ROL_MiniTensor_BoundConstraint_Def.hpp
#    function/ROL_MiniTensor_EqualityConstraint.hpp
#    function/ROL_MiniTensor_EqualityConstraint_Def.hpp
#    function/ROL_MiniTensor_Function.hpp
#    function/ROL_MiniTensor_Function_Def.hpp
#    function/ROL_MiniTensor_InequalityConstraint.hpp
#    function/ROL_MiniTensor_InequalityConstraint_Def.hpp
#  )
  TRIBITS_SET_AND_INC_DIRS( DIR ${CMAKE_CURRENT_SOURCE_DIR}/function )
  TRIBITS_INCLUDE_DIRECTORIES( ${DIR} )
  APPEND_GLOB(HEADERS ${DIR}/*.hpp)

ENDIF()

#
# solver headers
#

IF (ROL_ENABLE_Boost AND ROL_ENABLE_MiniTensor)
#  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sol)
#  APPEND_SET(HEADERS
#    sol/ROL_MiniTensor_MiniSolver.hpp
#    sol/ROL_MiniTensor_MiniSolver_Def.hpp
#  )

  TRIBITS_SET_AND_INC_DIRS( DIR ${CMAKE_CURRENT_SOURCE_DIR}/sol )
  TRIBITS_INCLUDE_DIRECTORIES( ${DIR} )
  APPEND_GLOB(HEADERS ${DIR}/*.hpp)

ENDIF()

INSTALL(FILES ${HEADERS} DESTINATION ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})
