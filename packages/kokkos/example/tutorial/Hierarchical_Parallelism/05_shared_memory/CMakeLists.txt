
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# This is a tutorial, not a test, so we don't ask CTest to run it.
TRIBITS_ADD_EXECUTABLE(
  tutorial_hierarchicalparallelism_02_shared_memory
  SOURCES shared_memory.cpp
  COMM serial mpi
  DEPLIBS kokkoscontainers kokkoscore
  )
