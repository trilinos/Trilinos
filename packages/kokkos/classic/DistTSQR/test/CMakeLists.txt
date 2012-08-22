INCLUDE(TribitsAddExecutableAndTest)

# Performance and accuracy test suite for TSQR::DistTsqr (which
# combines triangular factors from different MPI processes).
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TSQR_DistTsqr_Accuracy
  SOURCES Tsqr_TestDistTsqr.cpp
  COMM mpi
  ARGS "--verify --ncols=5 --explicit --implicit --real"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  DEPLIBS kokkos kokkosnodetsqr kokkosdisttsqr
)


INCLUDE(TribitsAddExecutableAndTest)

# Accuracy test for TSQR::Tsqr (the full TSQR implementation).
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TSQR_FullTsqr_Accuracy
  SOURCES Tsqr_TestFullTsqr.cpp
  COMM mpi
  ARGS "--numRowsLocal=100 --numCols=5 --testFactorExplicit --testReal"
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 4
  DEPLIBS kokkos kokkosnodetsqr kokkosdisttsqr
)
