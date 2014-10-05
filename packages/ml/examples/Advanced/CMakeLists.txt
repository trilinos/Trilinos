
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE(
  ml_analyze
  SOURCES ml_analyze.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  ml_filtering
  SOURCES ml_filtering.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  ml_user_aggr
  SOURCES ml_user_aggr.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  ml_with_ifpack
  SOURCES ml_with_ifpack.cpp
  COMM serial mpi
  )
