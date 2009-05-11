
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

FILE(GLOB HEADERS *.h)

FILE(GLOB SOURCES *.cpp)

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  gaasp
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )





