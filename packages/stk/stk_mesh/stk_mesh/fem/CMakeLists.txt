
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

# PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../)

#
# Core headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
EntityTypes.hpp
FieldDeclarations.hpp
FieldTraits.hpp
Stencils.hpp
TopologyDimensions.hpp
TopologyHelpers.hpp
  )

APPEND_SET(SOURCES
EntityTypes.cpp
FieldTraits.cpp
TopologyDimensions.cpp
TopologyHelpers.cpp
  )

IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  stk_mesh_fem
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
