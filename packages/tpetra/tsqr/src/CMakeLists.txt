# Set up this subpackage's automatically generated header file, that
# contains all the macros declaring the subpackage's CMake
# configuration options.
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(HEADERS "")
SET(SOURCES "")
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# Add this subpackage's generated header file to the list of header
# files to install.
APPEND_SET(HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h)

# Add headers and sources to the lists of this subpackage's headers
# resp. sources.
APPEND_GLOB(HEADERS ${DIR}/Tsqr*.hpp)
APPEND_GLOB(HEADERS ${DIR}/KokkosTSQR*.hpp)
APPEND_GLOB(SOURCES ${DIR}/Tsqr*.cpp)

TRIBITS_ADD_LIBRARY(
  kokkostsqr
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

#
# Make a trivial change to this comment if you add / remove a file either to
# / from this directory, or to / from the 'impl' subdirectory.  That ensures
# that running "make" will also rerun CMake in order to regenerate Makefiles.
#
# Behold: Another such change that I hath wrought, and another.
#
