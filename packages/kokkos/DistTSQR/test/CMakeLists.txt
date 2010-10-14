INCLUDE(PackageAddExecutableAndTest)

# Performance and accuracy test suite for TSQR::DistTsqr (which
# combines triangular factors from different MPI processes).
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  TSQR_DistTsqr_Accuracy
  SOURCES Tsqr_TestDistTsqr.cpp
  COMM mpi
  ARGS "--verify --ncols=5 --explicit --implicit --real"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  DEPLIBS kokkos kokkosnodetsqr kokkosdisttsqr
)
