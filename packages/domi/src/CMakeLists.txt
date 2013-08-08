# -*- cmake -*-

#
# Include the TriBITS library macros
INCLUDE(TribitsLibraryMacros)

#
# Package-specific configuration options
IF(Domi_ENABLE_Epetra)
  SET(HAVE_EPETRA ON)
ENDIF()

IF(Domi_ENABLE_Tpetra)
  SET(HAVE_TPETRA ON)
ENDIF()

#
# Set the Domi_ORDINAL_TYPE
SET(${PACKAGE_NAME}_ORDINAL_TYPE ptrdiff_t
  CACHE STRING
  "Select the data type for Domi_Ordinal (Domi::Ordinal)."
  )

#
# Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# Define the header and source files (and directories)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

APPEND_SET(HEADERS
  Domi_ConfigDefs.hpp
  Domi_Version.hpp
  Domi_Utils.hpp
  Domi_Exceptions.hpp
  Domi_Slice.hpp
  Domi_MDIterator.hpp
  Domi_MDRevIterator.hpp
  Domi_MDArray.hpp
  Domi_MDArrayView.hpp
  Domi_MDArrayRCP.hpp
  Domi_MDComm.hpp
  )

APPEND_SET(SOURCES
  Domi_Version.cpp
  Domi_Utils.cpp
  Domi_Slice.cpp
  Domi_MDComm.cpp
  )

#
# Define the targets for package's library(s)
TRIBITS_ADD_LIBRARY(
  domi
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
