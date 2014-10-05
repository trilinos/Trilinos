
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE(
  ml_example1d
  SOURCES ml_example1d.c place_holder.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  ml_example1dGS
  SOURCES ml_example1dGS.c place_holder.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  mlguide
  SOURCES mlguide.c place_holder.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  mlguide_par
  SOURCES mlguide_par.c place_holder.cpp
  COMM serial mpi
  )
