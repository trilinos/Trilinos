
INCLUDE(TribitsAddExecutableAndTest)


TRIBITS_ADD_TEST(
  UserOPUserMat
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "Solver converged"
  )
