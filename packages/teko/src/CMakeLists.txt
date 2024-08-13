

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_Config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_Config.h
  )

#
# Core files
#

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/NS)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Tpetra)

IF(TEKO_HAVE_EPETRA)
   TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Epetra)
ENDIF()

IF (ML_ENABLE_TekoSmoothers)
   TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/mlutils)
ENDIF()

APPEND_GLOB(HEADERS
  *.hpp
  NS/*.hpp
  Tpetra/*.hpp
  )

APPEND_GLOB(SOURCES
  *.cpp
  NS/*.cpp
  Tpetra/*.cpp
  )

IF(NOT TEKO_HAVE_EPETRA)
  LIST(REMOVE_ITEM SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/NS/Teko_ALOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Teko_DiagonalPreconditionerOp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Teko_ProbingPreconditionerFactory.cpp
  )
  LIST(REMOVE_ITEM HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/NS/Teko_ALOperator.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Teko_DiagonalPreconditionerOp.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Teko_ProbingPreconditionerFactory.hpp
  )
ELSE()
  APPEND_GLOB(HEADERS Epetra/*.hpp)
  APPEND_GLOB(SOURCES Epetra/*.cpp)
ENDIF()

IF (ML_ENABLE_TekoSmoothers)
   APPEND_GLOB(HEADERS mlutils/*.hpp)
   APPEND_GLOB(SOURCES mlutils/*.cpp)
ENDIF()

######################################

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  teko
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
