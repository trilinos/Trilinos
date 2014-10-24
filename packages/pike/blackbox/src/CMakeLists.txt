INCLUDE(TribitsLibraryMacros)
INCLUDE(TrilinosCreateClientTemplateHeaders)

#
# A) Package-specific configuration options
#

#TRIBITS_CONFIGURE_FILE(${PARENT_PACKAGE_NAME}_config.hpp)
#TRIBITS_CONFIGURE_FILE(${SUBPACKAGE_NAME}_config.hpp)
#TRIBITS_CONFIGURE_FILE(${SUBPACKAGE_FULLNAME}_config.hpp)
TRIBITS_CONFIGURE_FILE(Pike_BlackBox_config.hpp)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
IF (NOT ${PACKAGE_NAME}_ENABLE_EXPERIMENTAL)
  LIST(REMOVE_ITEM  HEADERS Pike_Solver_TransientStepper.hpp)
  LIST(REMOVE_ITEM  SOURCES Pike_Solver_TransientStepper.cpp)
ENDIF()
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  pike-blackbox
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
