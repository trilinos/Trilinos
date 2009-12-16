
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

FILE(GLOB HEADERS *.hpp *.h)

FILE(GLOB SOURCES *.cpp)

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  rythmos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )


#
# Add a test-only model library
#

IF (${PACKAGE_NAME}_ENABLE_TESTS)
  SET(HEADERS)
  SET(SOURCES)

  INCLUDE_DIRECTORIES( 
    ${PACKAGE_SOURCE_DIR}/test/SinCos
    ${PACKAGE_SOURCE_DIR}/test/PolynomialModel
    )

  APPEND_SET(HEADERS
    ${PACKAGE_SOURCE_DIR}/test/SinCos/SinCosModel.hpp
    ${PACKAGE_SOURCE_DIR}/test/PolynomialModel/PolynomialModel.hpp
    )

  APPEND_SET(SOURCES
    ${PACKAGE_SOURCE_DIR}/test/SinCos/SinCosModel.cpp
    ${PACKAGE_SOURCE_DIR}/test/PolynomialModel/PolynomialModel.cpp
    )

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Sacado)
  IF (${PACKAGE_NAME}_ENABLE_Sacado)
    INCLUDE_DIRECTORIES( 
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

    INCLUDE_DIRECTORIES( 
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
    INCLUDE_DIRECTORIES( 
      ${PACKAGE_SOURCE_DIR}/adapters/gaasp/example
      )
  
    APPEND_SET(HEADERS
      ${PACKAGE_SOURCE_DIR}/adapters/gaasp/example/LorenzModel.hpp
      )
  
    APPEND_SET(SOURCES
      ${PACKAGE_SOURCE_DIR}/adapters/gaasp/example/LorenzModel.cpp
      )

  ENDIF()  

  PACKAGE_ADD_LIBRARY(
    rythmos_test_models
    NOINSTALLHEADERS ${HEADERS}
    SOURCES ${SOURCES}
    TESTONLY
    DEPLIBS rythmos
    )

ENDIF()  



