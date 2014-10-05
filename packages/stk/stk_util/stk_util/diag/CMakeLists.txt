
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

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# Core headers
#

FILE(GLOB HEADERS *.hpp *.h)
#removing references to mpih which is not available in Trilinos
LIST(REMOVE_ITEM HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/mpih.hpp")

FILE(GLOB SOURCES *.cpp)
#removing references to mpih which is not available in Trilinos
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/mpih.cpp")
LIST(REMOVE_ITEM HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/Exception.hpp")

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  stk_util_diag
  DEPLIBS stk_util_parallel stk_util_env
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )


INSTALL(FILES ${HEADERS} DESTINATION
    ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_util/diag/)

