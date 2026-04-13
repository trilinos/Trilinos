
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
  SolverFactoryMaxIters
  SOURCES SolverFactoryMaxIters.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS 
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SolverFactoryMaxRestarts
  SOURCES SolverFactoryMaxRestarts.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS 
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SolverFactoryLossOfAccuracy
  SOURCES SolverFactoryLossOfAccuracy.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS 
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SolverFactoryBreakdownDetected
  SOURCES SolverFactoryBreakdownDetected.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS 
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SolverFactoryOrthonormFailure
  SOURCES SolverFactoryOrthonormFailure.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
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
