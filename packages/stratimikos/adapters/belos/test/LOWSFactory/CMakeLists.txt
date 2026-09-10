
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_BINARY_DIR}/../../src)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Thyra_Belos_StatusTest_UnitTests
  SOURCES
    Thyra_Belos_StatusTest_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


TRIBITS_ADD_EXECUTABLE(
  test_single_belos_thyra_solver_driver
  SOURCES
    test_single_belos_thyra_solver_driver.cpp
    test_single_belos_thyra_solver.cpp
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  test_single_belos_thyra_solver_driver
  POSTFIX_AND_ARGS_0 FourByFour
    --max-iters=4 --no-use-preconditioner --show-all-tests  --matrix-file=FourByFour.mtx
  POSTFIX_AND_ARGS_1 nos1_np
    --no-use-preconditioner --show-all-tests --num-rhs=1 --block-size=1
    --gmres-krylov-length=300 --max-fwd-error=1e-13 --max-resid=1e-13 --max-iters=1000
    --max-solution-error=1e-6 --matrix-file=nos1.mtx
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack2)
IF (${PACKAGE_NAME}_ENABLE_Ifpack2)

  # The test with 8 rhs vectors fails.

TRIBITS_ADD_TEST(
  test_single_belos_thyra_solver_driver
  POSTFIX_AND_ARGS_0 nos5_kl190
    --show-all-tests --gmres-krylov-length=190 --max-solution-error=1e-2 --max-fwd-error=1e-13 --max-resid=1e-6
    --matrix-file=nos5.mtx
  POSTFIX_AND_ARGS_1 nos1
    --show-all-tests --num-rhs=1 --block-size=1 --gmres-krylov-length=300
    --max-fwd-error=1e-12 --max-resid=1e-12 --max-iters=1000 --max-solution-error=1e-6
    --matrix-file=nos1.mtx
  # POSTFIX_AND_ARGS_2 nos1_nrhs8
  #   --show-all-tests --num-rhs=8 --block-size=8 --gmres-krylov-length=300
  #   --max-fwd-error=1e-12 --max-resid=1e-12 --max-iters=1000 --max-solution-error=1e-6
  #   --matrix-file=nos1.mtx
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

ENDIF()

  TRIBITS_ADD_TEST(
    test_single_belos_thyra_solver_driver
    POSTFIX_AND_ARGS_0 nos1_bs8
      --no-use-preconditioner --show-all-tests --num-rhs=8 --block-size=8
      --gmres-krylov-length=300 --max-fwd-error=1e-13 --max-resid=1e-13
      --max-iters=1000 --max-solution-error=1e-6
      --matrix-file=nos1.mtx
    COMM serial
    )


TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyBelosTestFiles
  SOURCE_FILES
  FourByFour.mtx
  nos1.mtx
  nos5.mtx
  EXEDEPS
    test_single_belos_thyra_solver_driver
)
