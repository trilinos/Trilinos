
INCLUDE(TribitsLibraryMacros)

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


FILE(GLOB HEADERS *.hpp)
FILE(GLOB SOURCES *.cpp)

INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR})

#
# Core headers
#


#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  stk_topology
  NOINSTALLHEADERS ${HEADERS}
  DEPLIBS stk_util_env 
  SOURCES ${SOURCES}
  )

INSTALL(FILES ${HEADERS} DESTINATION
    ${CMAKE_INSTALL_PREFIX}/include/stk_topology/)
