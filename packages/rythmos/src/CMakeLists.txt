

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

FILE(GLOB HEADERS *.hpp *.h)

FILE(GLOB SOURCES *.cpp)

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  rythmos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )


#
# Add a test-only model library
#

IF (${PACKAGE_NAME}_ENABLE_TESTS OR ${PACKAGE_NAME}_ENABLE_EXAMPLES)
  SET(HEADERS)
  SET(SOURCES)

  TRIBITS_INCLUDE_DIRECTORIES(
    ${PACKAGE_SOURCE_DIR}/test/SinCos
    ${PACKAGE_SOURCE_DIR}/test/PolynomialModel
    ${PACKAGE_SOURCE_DIR}/test/LogTime
    )

  APPEND_SET(HEADERS
    ${PACKAGE_SOURCE_DIR}/test/SinCos/SinCosModel.hpp
    ${PACKAGE_SOURCE_DIR}/test/PolynomialModel/PolynomialModel.hpp
    ${PACKAGE_SOURCE_DIR}/test/LogTime/LogTimeModel.hpp
    )

  APPEND_SET(SOURCES
    ${PACKAGE_SOURCE_DIR}/test/SinCos/SinCosModel.cpp
    ${PACKAGE_SOURCE_DIR}/test/PolynomialModel/PolynomialModel.cpp
    ${PACKAGE_SOURCE_DIR}/test/LogTime/LogTimeModel.cpp
    )

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Sacado)
  IF (${PACKAGE_NAME}_ENABLE_Sacado)
    TRIBITS_INCLUDE_DIRECTORIES(
      ${PACKAGE_SOURCE_DIR}/test/VanderPol
      )

    APPEND_SET(HEADERS
      ${PACKAGE_SOURCE_DIR}/test/VanderPol/VanderPolModel.hpp
      )

    APPEND_SET(SOURCES
      ${PACKAGE_SOURCE_DIR}/test/VanderPol/VanderPolModel.cpp
      )
  ENDIF()

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
  IF (${PACKAGE_NAME}_ENABLE_EpetraExt)

    TRIBITS_INCLUDE_DIRECTORIES(
      REQUIRED_DURING_INSTALLATION_TESTING
      ${PACKAGE_SOURCE_DIR}/../epetraext/example/model_evaluator/DiagonalTransient
      )

    APPEND_SET(HEADERS
      ${PACKAGE_SOURCE_DIR}/../epetraext/example/model_evaluator/DiagonalTransient/EpetraExt_DiagonalTransientModel.hpp
      )

    APPEND_SET(SOURCES
      ${PACKAGE_SOURCE_DIR}/../epetraext/example/model_evaluator/DiagonalTransient/EpetraExt_DiagonalTransientModel.cpp
      )

  ENDIF()

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_GAASP)
  IF (${PACKAGE_NAME}_ENABLE_GAASP)
    TRIBITS_INCLUDE_DIRECTORIES(
      ${PACKAGE_SOURCE_DIR}/adapters/gaasp/example
      )

    APPEND_SET(HEADERS
      ${PACKAGE_SOURCE_DIR}/adapters/gaasp/example/LorenzModel.hpp
      )

    APPEND_SET(SOURCES
      ${PACKAGE_SOURCE_DIR}/adapters/gaasp/example/LorenzModel.cpp
      )

  ENDIF()

  TRIBITS_ADD_LIBRARY(
    rythmos_test_models
    NOINSTALLHEADERS ${HEADERS}
    SOURCES ${SOURCES}
    TESTONLY
    DEPLIBS rythmos
    )

ENDIF()



