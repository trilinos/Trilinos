
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_BINARY_DIR}/../../src)

TRIBITS_ADD_EXECUTABLE(
  test_single_amesos2_tpetra_solver_driver
  SOURCES
    test_single_amesos2_tpetra_solver_driver.cpp
    test_single_amesos2_tpetra_solver.cpp
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  test_single_amesos2_tpetra_solver_driver
  POSTFIX_AND_ARGS_0 KLU2
  --show-all-tests --solver-type=KLU2 --verbose --matrix-file=orsirr_2.mtx
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

IF(Amesos2_ENABLE_SuperLUDist)
  TRIBITS_ADD_TEST(
    test_single_amesos2_tpetra_solver_driver
    POSTFIX_AND_ARGS_0 SuperLU_DIST
    --show-all-tests --solver-type=SuperLU_DIST --verbose --matrix-file=orsirr_2.mtx
    COMM serial mpi
    NUM_MPI_PROCS 1
    )
ENDIF()

IF(Amesos2_ENABLE_CSS_MKL)
  TRIBITS_ADD_TEST(
    test_single_amesos2_tpetra_solver_driver
    POSTFIX_AND_ARGS_0 CSSMKL
    --show-all-tests --solver-type=CSSMKL --verbose --matrix-file=orsirr_2.mtx
    COMM serial mpi
    NUM_MPI_PROCS 4
    )
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyAmesos2TestFiles
  DEST_FILES
    orsirr_2.mtx
  SOURCE_DIR
    "${PACKAGE_SOURCE_DIR}/../amesos2/test/matrices"
  EXEDEPS
    test_single_amesos2_tpetra_solver_driver
  )
