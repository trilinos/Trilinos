
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tpetra_bl_gmres_multidiag
  SOURCES test_bl_gmres_multidiag.cpp BelosLinearMultiShiftProblem.hpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_MultipleSolves
  SOURCES MultipleSolves.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS
  COMM serial mpi
  )
