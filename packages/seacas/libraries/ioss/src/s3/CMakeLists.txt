# Copyright(C) 2025 National Technology & Engineering Solutions
# of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
# NTESS, the U.S. Government retains certain rights in this software.
#
# See packages/seacas/LICENSE for details

SET(HEADERS "")
SET(SOURCES "")

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.h)
APPEND_GLOB(SOURCES ${DIR}/*.C)

TRIBITS_INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}/../"
  "${CMAKE_CURRENT_BINARY_DIR}/../"
  )

list(APPEND HEADERS
  "${CMAKE_CURRENT_BINARY_DIR}/ios3_export.h")

TRIBITS_ADD_LIBRARY(
  Ios3
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS
  Ioss
  )

generate_export_header(Ios3)
target_include_directories(Ios3
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>")

IF (BUILD_TESTING)
  ENABLE_TESTING()
ENDIF()
