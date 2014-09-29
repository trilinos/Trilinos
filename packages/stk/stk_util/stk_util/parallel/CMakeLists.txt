

#
# A) Package-specific configuration options
#

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

#
# Core headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

FILE(GLOB HEADERS *.hpp *.h)
FILE(GLOB SOURCES *.cpp)

#APPEND_SET(HEADERS
#BroadcastArg.hpp
#DistributedIndex.hpp
#ParallelComm.hpp
#Parallel.hpp
#ParallelIndex.hpp
#ParallelInputStream.hpp
#ParallelReduce.hpp
#  )
#
#APPEND_SET(SOURCES
#BroadcastArg.cpp
#DistributedIndex.cpp
#ParallelComm.cpp
#Parallel.cpp
#ParallelIndex.cpp
#ParallelInputStream.cpp
#ParallelReduce.cpp
#  )
#
#IF (TPL_ENABLE_MPI)
#ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  stk_util_parallel
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

INSTALL(FILES ${HEADERS} DESTINATION
    ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_util/parallel/)
