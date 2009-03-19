
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

ASSERT_DEFINED(Teuchos_ENABLE_Boost)
SET(HAVE_RYTHMOS_TEUCHOS_BOOST ${Teuchos_ENABLE_Boost})

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

FILE(GLOB HEADERS *.hpp *.h)

FILE(GLOB SOURCES *.cpp)

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

PACKAGE_ADD_LIBRARY(
  rythmos_test_models
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES}
  TESTONLY
  DEPLIBS rythmos
  )




