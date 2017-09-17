
INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# All headers
#

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  XROL.hpp
)

# utility
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/utility)
INCLUDE_DIRECTORIES(${DIR})
INCLUDE_DIRECTORIES(${DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

# vector headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/vector)
INCLUDE_DIRECTORIES(${DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

# function headers
#SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/function)
#INCLUDE_DIRECTORIES(${DIR}
#                    ${DIR}/objective)
#APPEND_GLOB(HEADERS ${DIR}/*.hpp
#                    ${DIR}/objective/*.hpp)
ADD_SUBDIRECTORY(function)


#IF (TPL_ENABLE_MPI)
#ENDIF()


#
# C) Define the targets for package's library(s)
# 

  TRIBITS_ADD_LIBRARY(
    xrol
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )
