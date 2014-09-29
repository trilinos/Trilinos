

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


FILE(GLOB HEADERS *.hpp *.tcc)
FILE(GLOB SOURCES *.cpp)

INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR})

#
# Core headers
#

FILE(GLOB DETAIL_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/topology_detail/*.hpp)

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
    ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_topology/)

INSTALL(FILES ${DETAIL_HEADERS} DESTINATION
    ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_topology/topology_detail)
