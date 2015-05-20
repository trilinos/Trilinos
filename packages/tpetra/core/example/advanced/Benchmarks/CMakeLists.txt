TRIBITS_ADD_EXECUTABLE(
  import
  SOURCES import.cpp
  COMM mpi
)

TRIBITS_ADD_EXECUTABLE(
  vector
  SOURCES vector.cpp
  COMM serial mpi
)
