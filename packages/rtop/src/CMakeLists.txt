
INCLUDE(TrilinosCreateClientTemplateHeaders)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_Config.h)

# Windows shared-library import/export for DATA symbols in rtop.dll
# (e.g. RTOpPack::transpMap). Consumers need __declspec(dllimport) or LNK2019.
SET(CURRENT_PACKAGE RTOP)
CONFIGURE_FILE("${Trilinos_SOURCE_DIR}/packages/Trilinos_DLLExportMacro.h.in"
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_DLLExportMacro.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

#
# src
#


TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

#
# src/interfaces
#

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/interfaces)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

#
# src/support
#

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/support)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
APPEND_GLOB(HEADERS ${DIR}/*.h)
APPEND_GLOB(SOURCES ${DIR}/*.c)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

#
# src/ops_lib
#

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/ops_lib)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

#
# src/lapack
#


TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/lapack)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

#
# Binary src
#

# Must glob the binary src dir last to get all of the generated headers
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_Config.h )
APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_DLLExportMacro.h)

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  rtop
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEFINES -DRTOP_LIB_EXPORTS_MODE
  )
