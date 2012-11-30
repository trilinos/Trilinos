
INCLUDE(TribitsLibraryMacros)

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.h)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

ASSERT_DEFINED(TeuchosCore_ENABLE_MPI)
IF (NOT TeuchosCore_ENABLE_MPI)

  LIST(REMOVE_ITEM  HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/Teuchos_DefaultMpiComm.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Teuchos_MpiReductionOpSetter.hpp
    )

  LIST(REMOVE_ITEM  SOURCES
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
