 
INCLUDE(TribitsLibraryMacros)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../src)

APPEND_GLOB(HEADERS *.hpp)
APPEND_GLOB(SOURCES *.cpp)

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  muelu-utils # should it be renamed muelu-common ? (it's not an extra utility lib, it's mandatory for the main lib)
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS muelu
  )

