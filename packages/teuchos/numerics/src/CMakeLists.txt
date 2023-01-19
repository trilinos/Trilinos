

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  teuchosnumerics
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEFINES -DTEUCHOSNUMERICS_LIB_EXPORTS_MODE
  )

#
# Make a trivial change to this comment if you add a file to this
# directory.  That ensures that running "make" will also rerun CMake
# in order to regenerate Makefiles.
#
