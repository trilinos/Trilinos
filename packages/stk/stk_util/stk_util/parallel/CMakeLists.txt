
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
BroadcastArg.hpp
DistributedIndex.hpp
ParallelComm.hpp
Parallel.hpp
ParallelIndex.hpp
ParallelInputStream.hpp
ParallelReduce.hpp
  )

APPEND_SET(SOURCES
BroadcastArg.cpp
DistributedIndex.cpp
ParallelComm.cpp
Parallel.cpp
ParallelIndex.cpp
ParallelInputStream.cpp
ParallelReduce.cpp
  )

IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  stk_util_parallel
  NOINSTALLHEADERS ${HEADERS}
  DEPLIBS stk_util_env
  SOURCES ${SOURCES}
  )

INSTALL(FILES ${HEADERS} DESTINATION
    ${CMAKE_INSTALL_PREFIX}/include/stk_util/parallel/)
