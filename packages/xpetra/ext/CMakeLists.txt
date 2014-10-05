# CMakeLists inspired by Tpetra/src/CMakeLists.txt

INCLUDE(TribitsLibraryMacros)
INCLUDE(TrilinosCreateClientTemplateHeaders)

#
# B) Define the header and source files (and directories)
#

#
# src
#

# Note about the use of wildcard in CMakeLists.txt: CMake dont't know
# when new files is added. You need to re-run CMake manually to
# regenerate a build system incorporating the new file.
# YOU MUST ALSO TOUCH A CMAKE CONFIGURATION FILE WHEN YOU PUSH THE NEW
# FILE TO FORCE THE RECONFIGURE ON OTHER PEOPLE'S BUILDS.

SET(HEADERS "")
APPEND_GLOB(HEADERS *.hpp)

SET(SOURCES mac.cpp) # Avoid empty library (MacOs 'ar' utility doesn't work on empty library).
APPEND_GLOB(SOURCES *.cpp)

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})

TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

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
