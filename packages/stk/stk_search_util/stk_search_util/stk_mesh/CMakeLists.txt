
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

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../)
#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../)
#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../stk_util)

#
# Core headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
    CreateBoundingBox.hpp
    PrintBoundingBox.hpp
    PrintEntityProc.hpp
  )

APPEND_SET(SOURCES
    CreateBoundingBox.cpp
    PrintEntityProc.cpp
  )

IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  stk_search_util
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

#INSTALL(FILES ${HEADERS} DESTINATION
#    ${CMAKE_INSTALL_PREFIX}/include/stk_search)
