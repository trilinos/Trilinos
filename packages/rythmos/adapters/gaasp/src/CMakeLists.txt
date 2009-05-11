
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_config.h )

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_GAASP_INTERNAL_TPL)
IF(${PACKAGE_NAME}_ENABLE_GAASP_INTERNAL_TPL)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../tpl_src)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  rythmosgaasp
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS rythmos gaasp
  )





