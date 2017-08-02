

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  directoryTest_Original
  SOURCES directoryTest_Original.cpp
  NUM_MPI_PROCS 4
  COMM mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  directoryTest_Relic
  SOURCES directoryTest_Relic.cpp
  NUM_MPI_PROCS 4
  COMM mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  directoryTest_Tpetra
  SOURCES directoryTest_Tpetra.cpp
  NUM_MPI_PROCS 4
  COMM mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  directoryTest_Kokkos
  SOURCES directoryTest_Kokkos.cpp
  NUM_MPI_PROCS 4
  COMM mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  directoryTest_KokkosSimple
  SOURCES directoryTest_KokkosSimple.cpp
  NUM_MPI_PROCS 4
  COMM mpi
)