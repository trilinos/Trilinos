# It's not necessary to run most of the tests below in an MPI build
# ("COMM mpi"), since only two of them (DistTsqr and FullTsqr) need to
# run on more than one MPI process.  However, it's useful to have the
# tests around in an MPI build, so we also build the tests there.  In
# an MPI build, only Process 0 in MPI_COMM_WORLD runs the tests; the
# other ranks are quieted.

ASSERT_DEFINED(TPL_ENABLE_CUDA)
ASSERT_DEFINED(Kokkos_ENABLE_Cuda)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_CUSOLVER)

IF (TPL_ENABLE_CUDA AND Kokkos_ENABLE_Cuda AND ${PACKAGE_NAME}_ENABLE_CUBLAS AND ${PACKAGE_NAME}_ENABLE_CUSOLVER)
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

# Performance and accuracy test suite for TSQR::KokkosNodeTsqr
# TRIBITS_ADD_TEST(
#   NodeTsqr
#   NAME KokkosNodeTsqr
#   COMM serial mpi
#   ARGS "--NodeTsqr=KokkosNodeTsqr --numRows=100000 --numCols=10"
#   STANDARD_PASS_OUTPUT
#   NUM_MPI_PROCS 1
#   )

# mfh 22 Dec 2014: Disable this test, since KokkosNodeTsqr no longer
# works with the new Kokkos Node types.
#
# Performance and accuracy test suite for TSQR::KokkosNodeTsqr
# ("generic" intranode parallel TSQR).  We pick an odd number of
# partitions to ensure correct results in that case, not just for
# powers of two (which everybody tests first).  The number of
# partitions is the maximum parallelism available in the algorithm,
# but it's up to the Kokkos Node implementation to decide what
# hardware resources to use (e.g., how many CPU cores, how many
# threads, ...).
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  KokkosNodeTsqr
#  SOURCES Tsqr_TestKokkosNodeTsqr.cpp
#  COMM serial mpi
#  ARGS "--verify --numRows=100000 --numCols=10 --numPartitions=7 --cacheSizeHint=50000 --contiguousCacheBlocks"
#  STANDARD_PASS_OUTPUT
#  NUM_MPI_PROCS 1
#  )

#
# Tests for the distributed-memory (MPI) part of TSQR.
#

# Performance and accuracy test suite for TSQR::DistTsqr (which
# combines triangular factors from different MPI processes).
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DistTsqr_Accuracy
  SOURCES Tsqr_TestDistTsqr.cpp
  COMM mpi
  ARGS "--verify --ncols=5 --explicit --implicit --real"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
)

# Accuracy test for TSQR::Tsqr (the full TSQR implementation).
TRIBITS_ADD_EXECUTABLE(
  FullTsqr
  SOURCES Tsqr_TestFullTsqr.cpp
  COMM mpi
  )

SET(TSQR_FULL_BASE_ARGS "--testFactorExplicit")
SET(TSQR_FULL_COMPLEX_BROKEN OFF)
IF(TSQR_FULL_COMPLEX_BROKEN)
  SET(TSQR_FULL_BASE_ARGS "${TSQR_FULL_BASE_ARGS} --noTestComplex")
ENDIF()

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
