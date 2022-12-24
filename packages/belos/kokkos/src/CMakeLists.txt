

#
# A) Package-specific configuration options
#

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

#
# Core Files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
   BelosKokkosAdapter.hpp
  )

#
# C) Define the targets for package's library(s)
#
#TRIBITS_ADD_LIBRARY(
#  beloskokkos
#  HEADERS ${HEADERS}
#  SOURCES ${SOURCES}
#  DEPLIBS belos
#  )

