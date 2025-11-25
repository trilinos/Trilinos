
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

ADD_SUBDIRECTORIES(
# Platform
  Map
  MapExtractor
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
  IO
  Matrix
  MatrixMatrix
  MatrixUtils
  )

# uncomment this to ensure Xpetra header files include all necessary dependencies
#ADD_SUBDIRECTORY(compile_tests)

ADD_SUBDIRECTORIES(
  BlockedMultiVector
)
