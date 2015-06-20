
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

ADD_SUBDIRECTORIES(
# Platform
  Map
# Distributor
# Directory
# ImportExport
# ImportExport2
  MultiVector
# BasicPerfTest
# Utils
# CrsGraph
  CrsMatrix

# Tests specifically written for Xpetra:
  Matrix
  MatrixMatrix
  Cloner
  )

# uncomment this to ensure Xpetra header files include all necessary dependencies
#ADD_SUBDIRECTORY(compile_tests)

IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE)
  ADD_SUBDIRECTORIES(
    BlockedCrsMatrix
  )
ENDIF()
