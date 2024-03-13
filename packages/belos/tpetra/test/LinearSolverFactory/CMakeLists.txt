
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SolverFactory
  SOURCES SolverFactory.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS 
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SolverFactoryNaN
  SOURCES SolverFactoryNaN.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS 
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CustomSolverFactory
  SOURCES CustomSolverFactory.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS 
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_LinearSolverFactory
  SOURCES LinearSolverFactory.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS 
  COMM serial mpi
  )
