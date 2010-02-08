
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_Config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_Config.h )

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/interfaces/operator_vector/ana/fundamental)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(HEADERS ${DIR}/*.cpp)

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/interfaces/operator_vector/ana/extended)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/support/operator_vector/client_support)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/support/operator_vector/adapter_support)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/support/operator_vector/adapter_support/serial)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/interfaces/operator_solve/ana/fundamental)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/interfaces/operator_solve/ana/extended)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/support/operator_solve/client_support)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/interfaces/nonlinear/model_evaluator/ana/fundamental)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/support/nonlinear/model_evaluator/client_support)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/interfaces/nonlinear/solvers/ana/fundamental)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/support/nonlinear/solvers/client_support)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

SET_AND_INC_DIRS(DIR ${PACKAGE_SOURCE_DIR}/example/operator_vector)
APPEND_SET(HEADERS ${DIR}/sillyModifiedGramSchmidt.hpp)

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  thyra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
