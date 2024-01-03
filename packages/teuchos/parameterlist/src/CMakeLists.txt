
#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  teuchosparameterlist
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEFINES -DTEUCHOSPARAMETERLIST_LIB_EXPORTS_MODE
  )

if (WIN32)
	target_link_libraries (teuchosparameterlist PUBLIC ws2_32)
endif ()
