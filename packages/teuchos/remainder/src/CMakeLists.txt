

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)
APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)


#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  teuchosremainder
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEFINES -DTEUCHOS_LIB_EXPORTS_MODE
  )

#
# Make a trivial change to this comment if you add / remove a file to
# / from this directory.  That ensures that running "make" will also
# rerun CMake in order to regenerate Makefiles.
#
