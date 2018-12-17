
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

# FIXME (mfh 12 Dec 2011) The tests that use preconditioners will fail
# if Trilinos was built without IFPACK support.  See e.g., line 101 of
# test_single_belos_thyra_solver.cpp.  Fix this by making execution of
# the test conditional on whether IFPACK was enabled.  I would do that
# myself, but I'm not sure how to split up this combined test into a
# list of conditionally enabled tests.
TRIBITS_ADD_TEST(
  test_single_belos_thyra_solver_driver
  POSTFIX_AND_ARGS_0 FourByFour
    --max-iters=4 --no-use-preconditioner --show-all-tests  --matrix-file=FourByFour.mtx
  POSTFIX_AND_ARGS_1 nos5_kl190
    --show-all-tests --gmres-krylov-length=190 --max-solution-error=1e-2 --max-resid=1e-6
    --matrix-file=nos5.mtx
  POSTFIX_AND_ARGS_2 nos1_np
    --no-use-preconditioner --show-all-tests --num-rhs=1 --block-size=1
    --gmres-krylov-length=300 --max-fwd-error=1e-13 --max-resid=1e-13 --max-iters=1000
    --max-solution-error=1e-6 --matrix-file=nos1.mtx
  POSTFIX_AND_ARGS_3 nos1
    --show-all-tests --num-rhs=1 --block-size=1 --gmres-krylov-length=300
    --max-fwd-error=1e-13 --max-resid=1e-13 --max-iters=1000 --max-solution-error=1e-6
    --matrix-file=nos1.mtx
  POSTFIX_AND_ARGS_4 nos1_nrhs8
    --show-all-tests --num-rhs=8 --block-size=8 --gmres-krylov-length=300
    --max-fwd-error=1e-13 --max-resid=1e-13 --max-iters=1000 --max-solution-error=1e-6
    --matrix-file=nos1.mtx
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

ASSERT_DEFINED(CMAKE_HOST_SYSTEM_NAME)
IF (CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")

  # This test produces lots of errors when run on Windows with Visual
  # C++ 2008.  This might point to a problem with Belos but I don't
  # have time to track this down right now.
  
  TRIBITS_ADD_TEST(
    test_single_belos_thyra_solver_driver
    POSTFIX_AND_ARGS_0 nos1_bs8
      --no-use-preconditioner --show-all-tests --num-rhs=8 --block-size=8
      --gmres-krylov-length=300 --max-fwd-error=1e-13 --max-resid=1e-13
      --max-iters=1000 --max-solution-error=1e-6
      --matrix-file=nos1.mtx
    COMM serial
    )

ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyBelosTestFiles
  DEST_FILES
    bcsstk01.mtx                              
    bcsstk02.mtx                              
    bcsstk04.mtx                              
    nos1.mtx                                      
    nos5.mtx                                      
    FourByFour.mtx                          
    Diagonal.mtx                              
    MissingADiagonal.mtx                              
    KheadK.mtx                                  
    KheadSorted.mtx
    Khead.triS  
  SOURCE_PREFIX
    "In_"
  SOURCE_DIR
    "${PACKAGE_SOURCE_DIR}/../amesos/test/Test_Basic"
  EXEDEPS
    test_single_belos_thyra_solver_driver
  )
