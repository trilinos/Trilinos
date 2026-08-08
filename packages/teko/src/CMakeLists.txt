#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_Config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS} ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_Config.h)

#
# Core files
#

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/NS)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Tpetra)

IF(TEKO_HAVE_EPETRA)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Epetra)
ENDIF()

APPEND_GLOB(
  HEADERS
  *.hpp
  NS/*.hpp
  Tpetra/*.hpp
)

APPEND_GLOB(
  SOURCES
  *.cpp
  NS/*.cpp
  Tpetra/*.cpp
)

IF(TEKO_HAVE_EPETRA)
  APPEND_GLOB(HEADERS Epetra/*.hpp)
  APPEND_GLOB(SOURCES Epetra/*.cpp)
ENDIF()

######################################

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  teko
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
)
