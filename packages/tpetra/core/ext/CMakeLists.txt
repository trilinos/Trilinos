
INCLUDE(TrilinosCreateClientTemplateHeaders)

#
# A) Package-specific configuration options
#

# TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

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
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

# Must glob the binary dir last to get all of the auto-generated headers
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  tpetraext
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS tpetra
  )

# TRIBITS_ADD_TEST_DIRECTORIES(test)
