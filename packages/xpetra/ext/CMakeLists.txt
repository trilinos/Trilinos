# CMakeLists inspired by Tpetra/src/CMakeLists.txt

INCLUDE(TribitsLibraryMacros)
INCLUDE(TrilinosCreateClientTemplateHeaders)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "") #*.hpp
#APPEND_GLOB(HEADERS *.hpp)

SET(SOURCES mac.cpp XpetraExt_MatrixMatrix.hpp) # Avoid empty library (MacOs 'ar' utility doesn't work on empty library).
APPEND_GLOB(SOURCES *.cpp)

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})

TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/../sup/Matrix)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/../sup/StridedMap)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/../sup/BlockedCrsMatrix)

# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  xpetra-ext
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS xpetra
  )
