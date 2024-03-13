# It's not necessary to run most of the tests below in an MPI build
# ("COMM mpi"), since only two of them (DistTsqr and FullTsqr) need to
# run on more than one MPI process.  However, it's useful to have the
# tests around in an MPI build, so we also build the tests there.  In
# an MPI build, only Process 0 in MPI_COMM_WORLD runs the tests; the
# other ranks are quieted.

ASSERT_DEFINED(TPL_ENABLE_CUDA)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_CUSOLVER)

IF (TPL_ENABLE_CUDA AND Kokkos_ENABLE_CUDA AND ${PACKAGE_NAME}_ENABLE_CUBLAS AND ${PACKAGE_NAME}_ENABLE_CUSOLVER)
  SET (TpetraTSQR_ENABLE_CUDA_TESTS ON)
ELSE ()
  SET (TpetraTSQR_ENABLE_CUDA_TESTS OFF)
ENDIF ()

IF (TpetraTSQR_ENABLE_CUDA_TESTS)
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CuSolver
  SOURCES CuSolver.cpp
  COMM serial mpi
  ARGS ""
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
)
ENDIF (TpetraTSQR_ENABLE_CUDA_TESTS)

# Performance and accuracy test suite for TSQR::Combine (which factors
# cache blocks and combines triangular factors).

TRIBITS_ADD_EXECUTABLE(
  Combine
  SOURCES Tsqr_TestCombine.cpp
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  Combine
  NAME Combine_100rows_5cols
  COMM serial mpi
  ARGS "--verify --numRows=100 --numCols=5"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  Combine
  NAME Combine_100rows_50cols
  COMM serial mpi
  ARGS "--verify --numRows=100 --numCols=50"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  Combine
  NAME Combine_10000rows_11cols
  COMM serial mpi
  ARGS "--verify --numRows=10000 --numCols=11"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

# This executable can test any NodeTsqr subclass that
# TSQR::NodeTsqrFactory can create.  It can check accuracy (--verify)
# and/or timing (--benchmark).  For both of these, it can compare with
# LAPACK.  Thus, this can serve as a check for your LAPACK
# implementation as well.  Run the executable with --help to see all
# the options.  It builds with or without MPI, but only runs with one
# MPI process.

TRIBITS_ADD_EXECUTABLE(
  NodeTsqr
  SOURCES Tsqr_TestNodeTsqr.cpp
  COMM serial mpi
  )

SET(TSQR_SEQUENTIALTSQR_COMPLEX_BROKEN ON)
SET(TSQR_SEQUENTIALTSQR_BASE_ARGS "--verify --NodeTsqr=SequentialTsqr")
IF(TSQR_SEQUENTIALTSQR_COMPLEX_BROKEN)
  SET(TSQR_SEQUENTIALTSQR_BASE_ARGS "${TSQR_SEQUENTIALTSQR_BASE_ARGS} --noTestComplex")
ELSE()
  SET(TSQR_SEQUENTIALTSQR_BASE_ARGS "${TSQR_SEQUENTIALTSQR_BASE_ARGS} --testComplex")
ENDIF()

TRIBITS_ADD_TEST(
  NodeTsqr
  NAME SequentialTsqr_contiguousCacheBlocks
  COMM serial mpi
  ARGS "${TSQR_SEQUENTIALTSQR_BASE_ARGS} --numRows=100000 --numCols=10 --cacheBlockSize=5000 --contiguousCacheBlocks"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  NodeTsqr
  NAME SequentialTsqr_noncontiguousCacheBlocks
  COMM serial mpi
  ARGS "${TSQR_SEQUENTIALTSQR_BASE_ARGS} --numRows=100000 --numCols=10 --cacheBlockSize=5000"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  NodeTsqr
  NAME CombineNodeTsqr
  COMM serial mpi
  ARGS "--verify --NodeTsqr=CombineNodeTsqr --numRows=1000 --numCols=15"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

IF (TpetraTSQR_ENABLE_CUDA_TESTS)
  TRIBITS_ADD_TEST(
    NodeTsqr
    NAME CuSolverNodeTsqr_11_5
    COMM serial mpi
    ARGS "--verify --NodeTsqr=CuSolverNodeTsqr --numRows=11 --numCols=5"
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_TEST(
    NodeTsqr
    NAME CuSolverNodeTsqr_5000_20
    COMM serial mpi
    ARGS "--verify --NodeTsqr=CuSolverNodeTsqr --numRows=5000 --numCols=20"
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1
    )
ENDIF ()

#
# Tests for the distributed-memory (MPI) part of TSQR.
#

# Performance and accuracy test suite for TSQR::DistTsqr (which
# combines triangular factors from different MPI processes).

# Accuracy test for TSQR::Tsqr (the full TSQR implementation).
TRIBITS_ADD_EXECUTABLE(
  DistTsqr
  SOURCES Tsqr_TestDistTsqr.cpp
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  DistTsqr
  NAME DistTsqr_1_proc
  COMM serial mpi
  ARGS "--verify --ncols=5 --explicit --implicit --real"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
)

TRIBITS_ADD_TEST(
  DistTsqr
  NAME DistTsqr_4_proc
  COMM mpi
  ARGS "--verify --ncols=5 --explicit --implicit --real"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 4
)

# Accuracy test for TSQR::Tsqr (the full TSQR implementation).
TRIBITS_ADD_EXECUTABLE(
  FullTsqr
  SOURCES Tsqr_TestFullTsqr.cpp
  COMM mpi
  )

SET(TSQR_FULL_BASE_ARGS "--testFactorExplicit")

TRIBITS_ADD_TEST(
  FullTsqr
  NAME FullTsqr_Accuracy_100rows_5cols
  COMM mpi
  ARGS "--numRowsLocal=100 --numCols=5 ${TSQR_FULL_BASE_ARGS}"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
  FullTsqr
  NAME FullTsqr_Accuracy_100rows_20cols
  COMM mpi
  ARGS "--numRowsLocal=100 --numCols=20 ${TSQR_FULL_BASE_ARGS}"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
  FullTsqr
  NAME FullTsqr_Accuracy_10000rows_5cols
  COMM mpi
  ARGS "--numRowsLocal=10000 --numCols=5 ${TSQR_FULL_BASE_ARGS}"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
  FullTsqr
  NAME FullTsqr_Accuracy_10000rows_20cols
  COMM mpi
  ARGS "--numRowsLocal=10000 --numCols=20 ${TSQR_FULL_BASE_ARGS}"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 4
)

IF (TpetraTSQR_ENABLE_CUDA_TESTS)
  TRIBITS_ADD_TEST(
    FullTsqr
    NAME FullTsqr_Accuracy_1000rows_15cols_CuSolver
    COMM mpi
    ARGS "--numRowsLocal=1000 --numCols=15 --NodeTsqr=CuSolverNodeTsqr ${TSQR_FULL_BASE_ARGS}"
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 4
  )
ENDIF ()

IF(TSQR_SEQUENTIALTSQR_COMPLEX_BROKEN)
  SET(TSQR_FULL_BASE_ARGS_SEQ "--noTestComplex")
ELSE()
  SET(TSQR_FULL_BASE_ARGS_SEQ "--testComplex")
ENDIF()

TRIBITS_ADD_TEST(
  FullTsqr
  NAME FullTsqr_Accuracy_5000rows_100cols_Sequential
  COMM mpi
  ARGS "--numRowsLocal=5000 --numCols=100 --NodeTsqr=SequentialTsqr ${TSQR_FULL_BASE_ARGS_SEQ}"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 4
)
