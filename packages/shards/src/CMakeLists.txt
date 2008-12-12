
INCLUDE(TrilinosPackageLibraryHelpers)

#
# A) Package-specific configuration options
#

TRILINOS_PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# CORE
APPEND_SET(HEADERS
  Shards_CellTopologyData.h
  Shards_Array.hpp
  Shards_ArrayVector.hpp
  Shards_CellTopology.hpp
  Shards_CellTopologyTraits.hpp
  Shards_BasicTopologies.hpp
  Shards_IndexList.hpp
  Shards_SimpleArrayOps.hpp
  Shards_CellTopology.hpp
  Shards_TypeList.hpp
  ) 

# CORE
APPEND_SET(SOURCES
  Shards_Array.cpp
  Shards_CellTopology.cpp
  Shards_BasicTopologies.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRILINOS_PACKAGE_ADD_LIBRARY(
  shards
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
