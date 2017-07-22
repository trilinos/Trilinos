

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  directoryTest_Original
  SOURCES directoryComparisonTest_Original.cpp
  NUM_MPI_PROCS 4
  COMM mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  directoryTest_Relic
  SOURCES directoryComparisonTest_Relic.cpp
  NUM_MPI_PROCS 4
  COMM mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  directoryTest_Tpetra
  SOURCES directoryComparisonTest_Tpetra.cpp
  NUM_MPI_PROCS 4
  COMM mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  directoryTest_Kokkos
  SOURCES directoryComparisonTest_Kokkos.cpp
  NUM_MPI_PROCS 4
  COMM mpi
)