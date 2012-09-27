INCLUDE(TribitsAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  Verbose
  SOURCES Verbose.cpp
  COMM serial mpi
  DEPLIBS muelu xpetra
  )
