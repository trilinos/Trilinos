

TRIBITS_ADD_EXECUTABLE(
  ml_maxwell
  SOURCES ml_maxwell.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  ml_simple_max
  SOURCES ml_simple_max.c place_holder.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  ml_read_maxwell
  SOURCES ml_read_maxwell.cpp
  COMM serial mpi
  )
