
INCLUDE(PackageLibraryMacros)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../src)

APPEND_GLOB(HEADERS *.hpp)
APPEND_GLOB(SOURCES *.cpp)

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  muelu-utils
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

