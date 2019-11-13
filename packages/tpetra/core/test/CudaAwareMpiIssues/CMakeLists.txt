TRIBITS_ADD_EXECUTABLE(
  randomMaps
  SOURCES randomMaps.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  repeatedMaps
  SOURCES repeatedMaps.cpp
  COMM serial mpi
  )
