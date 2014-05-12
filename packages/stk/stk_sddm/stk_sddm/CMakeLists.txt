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


#
# Core headers
#

INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR})
INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_BINARY_DIR}/stk_util/stk_util)

FILE(GLOB HEADERS *.hpp)

FILE(GLOB SOURCES *.cpp)

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  stk_sddm
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

