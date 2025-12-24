

#
# B) Define the header and source files (and include directories)
#
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(HEADERS "")
SET(SOURCES "")
LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h)

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.h)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

ASSERT_DEFINED(TeuchosCore_ENABLE_MPI)
IF (NOT TeuchosCore_ENABLE_MPI)

  LIST(REMOVE_ITEM  HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/Teuchos_Details_MpiCommRequest.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Teuchos_DefaultMpiComm.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Teuchos_MpiReductionOpSetter.hpp
    )

  LIST(REMOVE_ITEM  SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/Teuchos_Details_MpiCommRequest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Teuchos_DefaultMpiComm.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Teuchos_MpiReductionOpSetter.cpp
    )

ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  teuchoscomm
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEFINES -DTEUCHOSCOMM_LIB_EXPORTS_MODE
  )

#
# Make a trivial change to this comment if you add / remove a file
# either to / from this directory, or to / from the 'impl'
# subdirectory.  That ensures that running "make" will also rerun
# CMake in order to regenerate Makefiles.  Here is another such 
# change.
#
