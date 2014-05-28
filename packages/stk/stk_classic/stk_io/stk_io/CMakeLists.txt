

INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(STK_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

#
# Core headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../util/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/..)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../stk_util)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../seacas/ioss/src)

FILE(GLOB HEADERS *.hpp)
FILE(GLOB SOURCES *.cpp)

IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  stk_io
  NOINSTALLHEADERS ${HEADERS}
  DEPLIBS stk_mesh_base stk_mesh_fem
  SOURCES ${SOURCES}
  )

INSTALL(FILES ${HEADERS} DESTINATION
${CMAKE_INSTALL_PREFIX}/${Trilinos_INSTALL_INCLUDE_DIR}/stk_io/)

# AGS: I moved this to the bottom so that the link order was correct
# with  -lstk_io_uti before -lstk_io.  I think the DEPLIBS added in
# util/CMakeLists.txt should have done this, but it didn't work for me.
ADD_SUBDIRECTORY(util)
