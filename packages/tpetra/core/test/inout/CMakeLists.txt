TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixMarket_Tpetra_Issue11704
  SOURCES 
    Issue11704.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)


# This test should work for any number of MPI processes, as long as
# the sparse matrix to be read in has enough rows that it may be
# distributed over all of them, and as long as Proc 0 has enough
# memory to handle all the sparse matrix's data.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixMarket_Tpetra_CrsMatrix_InOutTest
  SOURCES
    MatrixMarket_Tpetra_CrsMatrix_InOutTest.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE(
  MatrixMarket_Tpetra_CrsMatrix_Dist
  SOURCES
    MatrixMarket_Tpetra_CrsMatrix_Dist.cpp
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_simple
  ARGS "--file=simple.mtx"
  COMM mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_simple
  ARGS "--file=simple.mtx"
  COMM mpi
  NUM_MPI_PROCS 3
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_simple
  ARGS "--file=simple.mtx"
  COMM mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_simple
  ARGS "--file=simple.mtx"
  COMM mpi
  NUM_MPI_PROCS 6
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_simple
  ARGS "--file=simple.mtx"
  COMM mpi
  NUM_MPI_PROCS 10
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_rmat
  ARGS "--file=rmat_11_lap.mtx"
  COMM mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_rmat
  ARGS "--file=rmat_11_lap.mtx"
  COMM mpi
  NUM_MPI_PROCS 3
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_rmat
  ARGS "--file=rmat_11_lap.mtx"
  COMM mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_rmat
  ARGS "--file=rmat_11_lap.mtx"
  COMM mpi
  NUM_MPI_PROCS 6
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_rmat
  ARGS "--file=rmat_11_lap.mtx"
  COMM mpi
  NUM_MPI_PROCS 10
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_rmat_nodiag
  ARGS "--file=rmat_11_lap.mtx --diagonal=exclude"
  COMM mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_rmat_nodiag
  ARGS "--file=rmat_11_lap.mtx --diagonal=exclude"
  COMM mpi
  NUM_MPI_PROCS 3
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_rmat_nodiag
  ARGS "--file=rmat_11_lap.mtx --diagonal=exclude"
  COMM mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_rmat_nodiag
  ARGS "--file=rmat_11_lap.mtx --diagonal=exclude"
  COMM mpi
  NUM_MPI_PROCS 6
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_rmat_nodiag
  ARGS "--file=rmat_11_lap.mtx --diagonal=exclude"
  COMM mpi
  NUM_MPI_PROCS 10
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE(
  MatrixMarket_Tpetra_CrsMatrix_Dist_Binary
  SOURCES
    MatrixMarket_Tpetra_CrsMatrix_Dist_Binary.cpp
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist_Binary
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_Binary_simple
  ARGS "--file=simple.mtx"
  COMM mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist_Binary
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_Binary_simple
  ARGS "--file=simple.mtx"
  COMM mpi
  NUM_MPI_PROCS 3
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist_Binary
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_Binary_simple
  ARGS "--file=simple.mtx"
  COMM mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist_Binary
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_Binary_simple
  ARGS "--file=simple.mtx"
  COMM mpi
  NUM_MPI_PROCS 6
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist_Binary
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_Binary_simple
  ARGS "--file=simple.mtx"
  COMM mpi
  NUM_MPI_PROCS 10
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist_Binary
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_Binary_rmat
  ARGS "--file=rmat_11_lap.mtx"
  COMM mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist_Binary
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_Binary_rmat
  ARGS "--file=rmat_11_lap.mtx"
  COMM mpi
  NUM_MPI_PROCS 3
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist_Binary
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_Binary_rmat
  ARGS "--file=rmat_11_lap.mtx"
  COMM mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist_Binary
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_Binary_rmat
  ARGS "--file=rmat_11_lap.mtx"
  COMM mpi
  NUM_MPI_PROCS 6
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist_Binary
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_Binary_rmat
  ARGS "--file=rmat_11_lap.mtx"
  COMM mpi
  NUM_MPI_PROCS 10
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist_Binary
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_BinaryPerProcess_simple
  ARGS "--file=simple.mtx --per-process"
  COMM mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist_Binary
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_BinaryPerProcess_simple
  ARGS "--file=simple.mtx --per-process"
  COMM mpi
  NUM_MPI_PROCS 3
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist_Binary
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_BinaryPerProcess_simple
  ARGS "--file=simple.mtx --per-process"
  COMM mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist_Binary
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_BinaryPerProcess_simple
  ARGS "--file=simple.mtx --per-process"
  COMM mpi
  NUM_MPI_PROCS 6
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist_Binary
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_BinaryPerProcess_simple
  ARGS "--file=simple.mtx --per-process"
  COMM mpi
  NUM_MPI_PROCS 10
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist_Binary
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_BinaryPerProcess_rmat
  ARGS "--file=rmat_11_lap.mtx --per-process"
  COMM mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist_Binary
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_BinaryPerProcess_rmat
  ARGS "--file=rmat_11_lap.mtx --per-process"
  COMM mpi
  NUM_MPI_PROCS 3
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist_Binary
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_BinaryPerProcess_rmat
  ARGS "--file=rmat_11_lap.mtx --per-process"
  COMM mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist_Binary
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_BinaryPerProcess_rmat
  ARGS "--file=rmat_11_lap.mtx --per-process"
  COMM mpi
  NUM_MPI_PROCS 6
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_Dist_Binary
  NAME MatrixMarket_Tpetra_CrsMatrix_Dist_BinaryPerProcess_rmat
  ARGS "--file=rmat_11_lap.mtx --per-process"
  COMM mpi
  NUM_MPI_PROCS 10
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixMarket_Tpetra_CrsGraph_InOutTest
  SOURCES
    MatrixMarket_Tpetra_CrsGraph_InOutTest.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )


  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MatrixMarket_Operator_Test
    SOURCES
      MatrixMarket_Tpetra_OperatorTest.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    ARGS ""
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

# This test should work for any number of MPI processes, as long as
# the dense matrix to be read in has enough rows that it may be
# distributed over all of them, and as long as Proc 0 has enough
# memory to handle all the dense matrix's data.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixMarket_Tpetra_MultiVector_InOutTest
  SOURCES
    MatrixMarket_Tpetra_MultiVector_InOutTest.cpp
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# This test has an upper bound on the number of MPI processes for
# which it works, since one of the test matrices has fixed dimensions.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixMarket_Tpetra_Map_InOutTest
  SOURCES
    MatrixMarket_Tpetra_Map_InOutTest.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
  )

# This test only makes sense for 1 MPI process.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Bug5800
  SOURCES
    Bug5800.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

# This test doesn't actually run as part of the check-in tests; it
# just builds.  It's good for testing Matrix Market input and output
# for Map and CrsMatrix, with arbitrary files.  The test should work
# with or without MPI support.
TRIBITS_ADD_EXECUTABLE(
  MatrixMarket_Tpetra_CrsMatrix_FileTest
  SOURCES
    MatrixMarket_Tpetra_CrsMatrix_FileTest.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  MatrixMarket_Tpetra_CrsGraph_FileTest
  SOURCES
    MatrixMarket_Tpetra_CrsGraph_FileTest.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  )

# This test has both lower and upper bounds on the number of MPI
# processes for which it makes sense to run.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Bug6288
  SOURCES Bug6288.cpp
  ARGS ""
  COMM mpi
  NUM_MPI_PROCS 2-32
  STANDARD_PASS_OUTPUT
  )

# This test only makes sense in an MPI build,
# with at least 2 MPI processes.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  PackTriples
  SOURCES
    PackTriples.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CooMatrixImpl
  SOURCES
    CooMatrixImpl.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# This test only makes sense in an MPI build,
# with at least 2 MPI processes.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CooMatrix
  SOURCES
    CooMatrix.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

# I want this test to _build_ in either MPI or non-MPI ("serial"),
# but I only want it to _run_ in an MPI build with exactly 2 MPI
# processes.
TRIBITS_ADD_EXECUTABLE(
  ReadTriples
  SOURCES
    ReadTriples.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  ReadTriples
  NAME ReadTriples_MPI_2
  ARGS ""
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

# Command that copies files to the executable directory.
TRIBITS_COPY_FILES_TO_BINARY_DIR(TpetraCore_inout_copy_files
  SOURCE_FILES
    simple.mtx
    rmat_11.mtx
    rmat_11_lap.mtx
    ReducedMatrix.mm
    ReducedMap.mm
)

TRIBITS_ADD_TEST(
  MatrixMarket_Tpetra_CrsMatrix_FileTest
  NAME ReadMapsWithGaps
  ARGS "--mapInputFilename=ReducedMap.mm --matrixInputFilename=ReducedMatrix.mm"
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
