INCLUDE(TrilinosCreateClientTemplateHeaders)

#
# A) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

IF(${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters)
  SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/thyra)
  APPEND_GLOB(HEADERS ${DIR}/*.hpp)
  APPEND_GLOB(SOURCES ${DIR}/*.cpp)
  TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})
ENDIF(${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters)

# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
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
