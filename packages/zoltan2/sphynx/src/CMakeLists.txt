
#
# A) Package-specific configuration options
#


TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)


#
# B) Define the header files, the source files, and the include directories
#

SET(HEADERS "")
SET(SOURCES "")


SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})

APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  ${PARAM_HEADER}
  )

# for Trilinos_version.h
ASSERT_DEFINED(PROJECT_BINARY_DIR)
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR})

# for Zoltan2_config.h
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  zoltan2sphynx
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
