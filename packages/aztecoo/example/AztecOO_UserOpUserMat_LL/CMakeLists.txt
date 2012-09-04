
INCLUDE(TribitsAddExecutableAndTest)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Poisson2d_LL
  SOURCES cxx_main.cpp Poisson2dOperator.cpp
  ARGS
    "21 21"
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "Solver converged"
  )
