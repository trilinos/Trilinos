INCLUDE(TrilinosCreateClientTemplateHeaders)

#
# A) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

IF(${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters)
  TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/thyra)
  APPEND_GLOB(HEADERS ${DIR}/*.hpp)
  APPEND_GLOB(SOURCES ${DIR}/*.cpp)
  TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})
ENDIF(${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters)

# Must glob the binary dir last to get all of the auto-generated headers
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#
# B) Define the targets for package's library(s)
#

IF(NOT "${SOURCES}" STREQUAL "")
  TRIBITS_ADD_LIBRARY(
    ifpack2-adapters
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS ifpack2
  )
ENDIF(NOT "${SOURCES}" STREQUAL "")

# If you add any .cpp files to this directory or its subdirectories
# that you want to add to the library, be sure to make a trivial
# change to this file.  This will ensure that CMake runs again.
# Here is another such change.
