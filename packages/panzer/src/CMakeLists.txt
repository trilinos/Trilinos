
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
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/dofmngr)

APPEND_SET(HEADERS
  Panzer_Version.hpp
  Panzer_ConfigDefs.hpp
  dofmngr/Panzer_FieldPattern.hpp
  dofmngr/Panzer_IntrepidFieldPattern.hpp
  dofmngr/Panzer_GeometricAggFieldPattern.hpp
  dofmngr/Panzer_FieldAggPattern.hpp
  dofmngr/Panzer_ConnManager.hpp
  ) 

APPEND_SET(SOURCES
  Panzer_Version.cpp
  dofmngr/Panzer_FieldPattern.cpp
  dofmngr/Panzer_IntrepidFieldPattern.cpp
  dofmngr/Panzer_GeometricAggFieldPattern.cpp
  dofmngr/Panzer_FieldAggPattern.cpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  panzer
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

