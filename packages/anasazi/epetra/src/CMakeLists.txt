
INCLUDE(TrilinosPackageLibraryHelpers)

#
# A) Package-specific configuration options
#

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

#
# Core Files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  AnasaziEpetraAdapter.hpp 
  )

APPEND_SET(SOURCES
  AnasaziEpetraAdapter.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRILINOS_PACKAGE_ADD_LIBRARY(
  anasaziepetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS anasazi
  )

#
# D) Export the dependancy variables of this package for other packages to use
#

TRILINOS_PACKAGE_EXPORT_DEPENDENCY_VARIABLES()
