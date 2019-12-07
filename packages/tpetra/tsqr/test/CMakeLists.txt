# It's not necessary to run the first five tests in an MPI build
# ("COMM mpi"), since none of them need to run on more than one MPI
# process.  However, it's useful to have the tests around in an MPI
# build, so we also build the tests there.  In an MPI build, only
# Process 0 in MPI_COMM_WORLD runs the tests; the other ranks are
# quieted.

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

# Test NodeTsqrFactory and NodeTsqr subclasses generically.
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

# Performance and accuracy test suite for TSQR::KokkosNodeTsqr
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  KokkosHostTsqr
  SOURCES Tsqr_TestKokkosNodeTsqr.cpp
  COMM serial mpi
  ARGS "--verify --numRows=100000 --numCols=10"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

# This test uses LAPACK's QR factorization to get a reference for
# performance and accuracy.  It doesn't run any parts of the TSQR
# algorithm, but it does depend on some TSQR test code (for generating
# the test matrix and measuring accuracy).
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Lapack
  SOURCES Tsqr_TestLapack.cpp
  COMM serial mpi
  ARGS "--verify --nrows=1000 --ncols=10 --ntrials=10"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

# Performance and accuracy test suite for TSQR::TBB::TbbTsqr
# (shared-memory parallel cache-blocked TSQR, parallelized via Intel's
# Threading Building Blocks library).
#
# Only build TBB-enabled TSQR if (surprise!) TBB is enabled.
IF (KokkosTSQR_ENABLE_TBB)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TbbTsqr
    SOURCES Tsqr_TestTbbTsqr.cpp
    COMM serial mpi
    ARGS "--verify --nrows=100000 --ncols=10 --cache-block-size=50000 --contiguous-cache-blocks"
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1
    )
ENDIF()

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

SET(TSQR_FULL_COMPLEX_BROKEN ON)
SET(TSQR_FULL_KOKKOSNODETSQR_BROKEN ON)
SET(TSQR_FULL_BASE_ARGS "--testFactorExplicit")
IF(TSQR_FULL_COMPLEX_BROKEN)
  SET(TSQR_FULL_BASE_ARGS "${TSQR_FULL_BASE_ARGS} --noTestComplex")
ENDIF()
IF(TSQR_FULL_KOKKOSNODETSQR_BROKEN)
  SET(TSQR_FULL_BASE_ARGS "${TSQR_FULL_BASE_ARGS} --NodeTsqr=SequentialTsqr")
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
