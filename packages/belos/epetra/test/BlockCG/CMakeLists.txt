

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  pseudo_cg_indefinite
  SOURCES test_pseudo_cg_indefinite.cpp   
  COMM serial mpi
  ARGS
    "--debug --verbose --tol=1e-5"
  STANDARD_PASS_OUTPUT
  )

