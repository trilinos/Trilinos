# -*- cmake -*-

#
# Include the TriBITS library macros
INCLUDE(TribitsLibraryMacros)

#
# Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(${PACKAGE_NAME}_ORDINAL_TYPE ptrdiff_t
  CACHE STRING
  "Select the data type for Domi_Ordinal (Domi::Ordinal)."
  )

#
# Define the header and source files (and directories)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

APPEND_SET(HEADERS
  Domi_MDArray.hpp
  Domi_MDArrayRCP.hpp
  Domi_MDArrayView.hpp
  Domi_Utils.hpp
  Domi_Slice.hpp
  )

APPEND_SET(SOURCES
  Domi_Version.cpp
  Domi_Utils.cpp
  Domi_Slice.cpp
  )

#
# Define the targets for package's library(s)
TRIBITS_ADD_LIBRARY(
  domi
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
