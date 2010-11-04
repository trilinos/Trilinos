
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

# PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

# SET(HEADERS ${HEADERS}
#   ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
#   )

#
# Core Files
#

# general headers, sources
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(HEADERS ${DIR}/*.cpp)
# PACKAGE_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})
# Must glob the binary dir last to get all of the auto-generated headers
# SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
# APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  tpetrarti
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS tpetra
  )

# PACKAGE_ADD_TEST_DIRECTORIES(test)
