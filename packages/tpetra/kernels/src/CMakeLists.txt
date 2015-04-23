
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------------

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/impl)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h)

#-----------------------------------------------------------------------------

TRIBITS_ADD_LIBRARY(
  tpetrakernels
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

#
# Make a trivial change here if you add a file either to this
# directory, or to the impl/ subdirectory.  Keep making changes
# whenever you add or remove a .cpp or .hpp file, and everything will
# be OK.  Oh look, it's another change.
#

