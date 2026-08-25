
INCLUDE(TrilinosCreateClientTemplateHeaders)

SET(HEADERS "")
SET(SOURCES "")

TRIBITS_SET_AND_INC_DIRS(srcDir ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${srcDir}/*.hpp)
APPEND_GLOB(SOURCES ${srcDir}/*.cpp)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${srcDir})

IF(${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters)
TRIBITS_SET_AND_INC_DIRS(tpetraSrcDir ${CMAKE_CURRENT_SOURCE_DIR}/tpetra)

APPEND_GLOB(HEADERS ${tpetraSrcDir}/*.hpp)
APPEND_GLOB(SOURCES ${tpetraSrcDir}/*.cpp)
ENDIF()

APPEND_GLOB(HEADERS ${CMAKE_CURRENT_BINARY_DIR}/*.hpp)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../../src)

# Add Belos solvers directory to include path for package-by-package builds
# This allows including Belos headers without the solvers/ prefix
IF(${PACKAGE_NAME}_ENABLE_Belos)
  # Check if this is a package build (Belos is external) or monolithic build
  IF(TPL_ENABLE_Belos)
    # Package-by-package build: Belos headers may be in install directory
    # The TPL include directories will be added automatically
  ELSE()
    # Monolithic build: add the source directory
    TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../../belos/tpetra/src/solvers)
  ENDIF()
ENDIF()

TRIBITS_ADD_LIBRARY(
  stratimikosbelos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
