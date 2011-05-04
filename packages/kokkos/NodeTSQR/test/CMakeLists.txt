INCLUDE(PackageAddExecutableAndTest)

# We do build MPI versions of the Combine, SequentialTsqr, LAPACK,
# and TbbTsqr tests.  However, they just quiet all processes but
# Rank 0 and do the tests there.

# Performance and accuracy test suite for TSQR::Combine (which
# factors cache blocks and combines triangular factors).
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  TSQR_Combine
  SOURCES Tsqr_TestCombine.cpp
  COMM serial mpi
  ARGS "--verify --nrows=100 --ncols=5"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  DEPLIBS kokkos kokkosnodetsqr
  )

# Performance and accuracy test suite for TSQR::SequentialTsqr
# (sequential cache-blocked TSQR).
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  TSQR_SequentialTsqr
  SOURCES Tsqr_TestSeqTsqr.cpp
  COMM serial mpi
  ARGS "--verify --nrows=100000 --ncols=10 --cache-block-size=50000 --contiguous-cache-blocks"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  DEPLIBS kokkos kokkosnodetsqr
  )

# This test uses LAPACK's QR factorization to get a reference for
# performance and accuracy.  It doesn't run any parts of the TSQR
# algorithm, but it does depend on some TSQR test code (for
# generating the test matrix and measuring accuracy).
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  TSQR_Lapack
  SOURCES Tsqr_TestLapack.cpp
  COMM serial mpi
  ARGS "--verify --nrows=1000 --ncols=10 --ntrials=10"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  DEPLIBS kokkos kokkosnodetsqr
  )

# Performance and accuracy test suite for TSQR::TBB::TbbTsqr
# (shared-memory parallel cache-blocked TSQR, parallelized via
# Intel's Threading Building Blocks library).
#
# Only build TBB-enabled TSQR if (surprise!) TBB is enabled.
IF (TSQR_ENABLE_Intel_TBB)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    TSQR_TbbTsqr
    SOURCES Tsqr_TestTbbTsqr.cpp
    COMM serial mpi
    ARGS "--verify --nrows=100000 --ncols=10 --cache-block-size=50000 --contiguous-cache-blocks"
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1
    DEPLIBS kokkos kokkosnodetsqr
    )
ENDIF()


# Performance and accuracy test suite for TSQR::KokkosNodeTsqr
# ("generic" intranode parallel TSQR).
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  TSQR_KokkosNodeTsqr
  SOURCES Tsqr_TestKokkosNodeTsqr.cpp
  COMM serial mpi
  ARGS "--verify --nrows=100000 --ncols=10 --cache-block-size=50000 --contiguous-cache-blocks"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  DEPLIBS kokkos kokkosnodetsqr
  )

