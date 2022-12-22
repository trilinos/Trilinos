


#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)


CONFIGURE_FILE("${Trilinos_SOURCE_DIR}/packages/Trilinos_DLLExportMacro.h.in"
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_DLLExportMacro.h)


#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")


# Set the base Trilinos directory in order to find Trilinos_version.h

ASSERT_DEFINED(PROJECT_SOURCE_DIR)
ASSERT_DEFINED(PROJECT_BINARY_DIR)

TRIBITS_INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})

#
# Include packages/src directory
#
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})


APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Set sources

APPEND_SET(SOURCES
  dummy.cpp
  )

TRIBITS_ADD_LIBRARY(
  shylu_ddcommon
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
