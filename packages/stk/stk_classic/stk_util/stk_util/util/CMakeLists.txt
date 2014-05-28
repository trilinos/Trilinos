
INCLUDE(TribitsLibraryMacros)

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

FILE(GLOB HEADERS *.hpp *.h)
FILE(GLOB SOURCES *.cpp)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS}/../)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS}/detail/)

#
# Core headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  stk_util_util
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

INSTALL(FILES ${HEADERS} DESTINATION
    ${CMAKE_INSTALL_PREFIX}/${Trilinos_INSTALL_INCLUDE_DIR}/stk_util/util/)
