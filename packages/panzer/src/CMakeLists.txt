
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.hpp
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# FAD
APPEND_SET(HEADERS
  Panzer_Version.hpp
  Panzer_ConfigDefs.hpp
  dofmngr_v2/Panzer_FieldPattern.hpp
  dofmngr_v2/Panzer_IntrepidFieldPattern.hpp
  dofmngr_v2/Panzer_GeometricAggFieldPattern.hpp
  dofmngr_v2/Panzer_FieldAggPattern.hpp
  ) 

APPEND_SET(SOURCES
  Panzer_Version.cpp
  dofmngr_v2/Panzer_FieldPattern.cpp
  dofmngr_v2/Panzer_IntrepidFieldPattern.cpp
  dofmngr_v2/Panzer_GeometricAggFieldPattern.cpp
  dofmngr_v2/Panzer_FieldAggPattern.cpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  panzer
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

