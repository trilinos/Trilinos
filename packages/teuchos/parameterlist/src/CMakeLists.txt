

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

IF (NOT TeuchosCore_ENABLE_yaml-cpp)

  LIST(REMOVE_ITEM  HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/Teuchos_YamlParameterListCoreHelpers.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Teuchos_YamlParser_decl.hpp
    )

 LIST(REMOVE_ITEM  SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/Teuchos_YamlParser.cpp
    )
ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  teuchosparameterlist
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEFINES -DTEUCHOSPARAMETERLIST_LIB_EXPORTS_MODE
  )
