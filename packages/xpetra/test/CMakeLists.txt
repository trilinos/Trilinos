INCLUDE(AddSubdirectories)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

ADD_SUBDIRECTORIES(
# Platform
  Map
# Distributor
# Directory
# ImportExport
# ImportExport2
#  MultiVector   ##### temporarily disabling this test
# BasicPerfTest
# Utils
# CrsGraph
  CrsMatrix
# BlockCrsGraph
# BlockMap
# BlockMultiVector
# VbrMatrix

# Tests specifically written for Xpetra:
  Operator
  )

# uncomment this to ensure Xpetra header files include all necessary dependencies
#ADD_SUBDIRECTORY(compile_tests)

IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE)
  ADD_SUBDIRECTORIES(
    BlockedCrsOperator
  )
ENDIF()
