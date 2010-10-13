# CMakeLists inspired by Tpetra/src/CMakeLists.txt

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

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})

# JG TODO: Is it the good macro to call ?
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Def)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Map)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Comm)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/VbrMatrix)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/MultiVector)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/BlockMap)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Platform)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/CrsMatrix)

APPEND_GLOB(HEADERS ${DIR}/*/*.cpp)
APPEND_GLOB(HEADERS ${DIR}/*/*.hpp)

PACKAGE_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
# APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_config.h )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  cthulhu
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )