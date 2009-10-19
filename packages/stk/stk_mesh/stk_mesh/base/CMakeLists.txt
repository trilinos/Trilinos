
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
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../stk_util)

#
# Core headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
Bucket.hpp
BulkData.hpp
Comm.hpp
DataTraitsClass.hpp
DataTraitsEnum.hpp
DataTraits.hpp
EntityComm.hpp
Entity.hpp
EntityKey.hpp
FieldData.hpp
Field.hpp
FieldParallel.hpp
GetBuckets.hpp
GetEntities.hpp
Ghosting.hpp
MetaData.hpp
Part.hpp
Property.hpp
Relation.hpp
Selector.hpp
SetOwners.hpp
Types.hpp
  )

APPEND_SET(SOURCES
Bucket.cpp
BulkData.cpp
BulkDataEndSync.cpp
BulkDataGhosting.cpp
BulkDataOwner.cpp
BulkDataParallelVerify.cpp
BulkDataSharing.cpp
Comm.cpp
DataTraits.cpp
EntityComm.cpp
Entity.cpp
Field.cpp
FieldData.cpp
FieldParallel.cpp
GetBuckets.cpp
GetEntities.cpp
MetaData.cpp
Part.cpp
Property.cpp
Relation.cpp
Selector.cpp
  )

IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  stk_mesh_base
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
