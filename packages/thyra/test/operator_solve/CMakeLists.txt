INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE(
  test_linear_op_with_solve
  SOURCES test_linear_op_with_solve.cpp
  COMM serial mpi
  )


IF (WIN32 AND NOT CYGWIN)
  SET(EPS_SCALE 2e+3)
ELSE()
  SET(EPS_SCALE 2e+2)
ENDIF()

PACKAGE_ADD_TEST(
  test_linear_op_with_solve
  NAME test_linear_op_with_solve
  POSTFIX_AND_ARGS_0 n1
    --n=1 --eps-scale=${EPS_SCALE}
  POSTFIX_AND_ARGS_1 n2
    --n=2 --eps-scale=${EPS_SCALE}
  NUM_MPI_PROCS 1
  COMM serial mpi
  )


IF (WIN32 AND NOT CYGWIN)
  SET(EPS_SCALE 2e+4)
ELSE()
  SET(EPS_SCALE 2e+2)
ENDIF()

PACKAGE_ADD_TEST(
  test_linear_op_with_solve
  NAME test_linear_op_with_solve
  POSTFIX_AND_ARGS_0 n4
    --n=4 --eps-scale=${EPS_SCALE}
  NUM_MPI_PROCS 1
  COMM serial mpi
  XHOSTTYPE Darwin 
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  DelayedLinearOpWithSolveUnitTests
  SOURCES
    DelayedLinearOpWithSolve_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  XHOSTTYPE Darwin 
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  DefaultMultiVectorLinearOpWithSolveUnitTests
  SOURCES
    DefaultMultiVectorLinearOpWithSolve_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  XHOSTTYPE Darwin 
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  DefaultBlockedTriangularLinearOpWithSolveUnitTests
  SOURCES
    DefaultBlockedTriangularLinearOpWithSolve_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  XHOSTTYPE Darwin 
  )


# Bug: 4699: ToDo: Fix the randomize(...) function to create random columns.
# Or, you can just generate a matrix with 1, 2, 3, 4 ... which should be a
# nicely conditioned matrix for small matrices.
