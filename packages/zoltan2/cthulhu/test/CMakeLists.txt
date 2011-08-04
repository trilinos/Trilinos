INCLUDE(AddSubdirectories)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

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
# CrsMatrix
# BlockCrsGraph
  BlockMap
# BlockMultiVector
# VbrMatrix

# Tests specifically written for Cthulhu:
  Operator
  )
