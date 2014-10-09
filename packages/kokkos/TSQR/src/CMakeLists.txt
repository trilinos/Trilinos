# Set up this subpackage's automatically generated header file, that
# contains all the macros declaring the subpackage's CMake
# configuration options.
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(HEADERS "")
SET(SOURCES "")
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# Add this subpackage's generated header file to the list of header
# files to install.
APPEND_SET(HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h)

# If TBB (Intel's Threading Building Blocks) is enabled, add headers
# and sources for TBB-enabled shared-memory parallel TSQR to the
# lists of this subpackage's headers resp. sources.
IF (${PACKAGE_NAME}_ENABLE_TBB)
  APPEND_GLOB(HEADERS ${DIR}/TbbTsqr*.hpp)
  APPEND_GLOB(SOURCES ${DIR}/TbbTsqr*.cpp)
ENDIF ()

# Add all other headers and sources (those not related to TBB) to the
# lists of this subpackage's headers resp. sources.
APPEND_GLOB(HEADERS ${DIR}/Tsqr*.hpp)
APPEND_GLOB(SOURCES ${DIR}/Tsqr*.cpp)

TRIBITS_ADD_LIBRARY(
  kokkostsqr
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

