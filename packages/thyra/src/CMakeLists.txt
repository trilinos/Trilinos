
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

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_Config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS *.hpp)
APPEND_GLOB(SOURCES *.cpp)

#
# src/interfaces/operator_vector/ana/fundamental
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/interfaces/operator_vector/ana/fundamental)
APPEND_GLOB(HEADERS interfaces/operator_vector/ana/fundamental/*.hpp)
APPEND_GLOB(HEADERS interfaces/operator_vector/ana/fundamental/*.cpp)

#
# src/interfaces/operator_vector/ana/extended
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/interfaces/operator_vector/ana/extended)
APPEND_GLOB(HEADERS interfaces/operator_vector/ana/extended/*.hpp)

#
# src/support/operator_vector/client_support
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/support/operator_vector/client_support)
APPEND_GLOB(HEADERS support/operator_vector/client_support/*.hpp)
APPEND_GLOB(SOURCES support/operator_vector/client_support/*.cpp)

#
# src/support/operator_vector/adapter_support
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/support/operator_vector/adapter_support)
APPEND_GLOB(HEADERS support/operator_vector/adapter_support/*.hpp)
APPEND_GLOB(SOURCES support/operator_vector/adapter_support/*.cpp)

#
# src/support/operator_vector/adapter_support/serial
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/support/operator_vector/adapter_support/serial)
APPEND_GLOB(HEADERS support/operator_vector/adapter_support/serial/*.hpp)

#
# src/interfaces/operator_solve/ana/fundamental
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/interfaces/operator_solve/ana/fundamental)
APPEND_GLOB(HEADERS interfaces/operator_solve/ana/fundamental/*.hpp)

#
# src/interfaces/operator_solve/ana/extended
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/interfaces/operator_solve/ana/extended)
APPEND_GLOB(HEADERS interfaces/operator_solve/ana/extended/*.hpp)

#
# src/support/operator_solve/client_support
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/support/operator_solve/client_support)
APPEND_GLOB(HEADERS support/operator_solve/client_support/*.hpp)
APPEND_GLOB(SOURCES support/operator_solve/client_support/*.cpp)

#
# src/interfaces/nonlinear/model_evaluator/ana/fundamental
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/interfaces/nonlinear/model_evaluator/ana/fundamental)
APPEND_GLOB(HEADERS interfaces/nonlinear/model_evaluator/ana/fundamental/*.hpp)
APPEND_GLOB(SOURCES interfaces/nonlinear/model_evaluator/ana/fundamental/*.cpp)

#
# src/support/nonlinear/model_evaluator/client_support
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/support/nonlinear/model_evaluator/client_support)
APPEND_GLOB(HEADERS support/nonlinear/model_evaluator/client_support/*.hpp)

#
# src/interfaces/nonlinear/solvers/ana/fundamental
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/interfaces/nonlinear/solvers/ana/fundamental)
APPEND_GLOB(HEADERS interfaces/nonlinear/solvers/ana/fundamental/*.hpp)

#
# src/support/nonlinear/solvers/client_support
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/support/nonlinear/solvers/client_support)
APPEND_GLOB(HEADERS support/nonlinear/solvers/client_support/*.hpp)

#
# thyra/example/operator_vector
#

INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/example/operator_vector)
APPEND_SET(HEADERS
  ${PACKAGE_SOURCE_DIR}/example/operator_vector/sillyModifiedGramSchmidt.hpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  thyra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
