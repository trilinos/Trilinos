INCLUDE(TribitsLibraryMacros)

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# Create a TriKota Version library first to trigger the regular TriBITS
# library initializations
#

TRIBITS_ADD_LIBRARY(
  trikotaversion
  HEADERS
    ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
    TriKota_ConfigDefs.hpp
    TriKota_Version.hpp
  SOURCES
    TriKota_Version.cpp
  )

#
# Build the TriKota library
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${Dakota_INCLUDE_DIRS}
  )

APPEND_SET(HEADERS
    TriKota_DirectApplicInterface.hpp
    TriKota_ThyraDirectApplicInterface.hpp
    TriKota_Driver.hpp
  ) 

APPEND_SET(SOURCES
    TriKota_DirectApplicInterface.cpp
    TriKota_ThyraDirectApplicInterface.cpp
    TriKota_Driver.cpp
  )

# Assume Dakota built with CMake

# Needed to find Dakota components built with ExternalProject_Add
# and their dependencies; only defined in recent versions of Dakota
IF(DEFINED Dakota_LINK_DIRS)
  LINK_DIRECTORIES(${Dakota_LINK_DIRS})
ENDIF()

# Do what TRIBITS_ADD_LIBRARY() would have done
PRINT_VAR(${PACKAGE_NAME}_LIBRARIES)
PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIBRARIES ${Dakota_LIBRARIES})
PRINT_VAR(${PACKAGE_NAME}_LIBRARIES)

# The DEPLIBS part is to tell CMake that trikota depends on the specified
# Dakota libraries.  The IMPORTEDLIBS part is to get the TRIBITS_ADD_LIBRARY
# macro to put the Dakota libraries TriKota_LIBRARIES when writing
# TriKotaConfig.cmake and Makefile.export.TriKota.  There should be a better
# way of doing this...
TRIBITS_ADD_LIBRARY(
  trikota
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS trikotaversion ${Dakota_LIBRARIES}
  IMPORTEDLIBS  ${Dakota_EXTRA_LIBRARIES}
  )
