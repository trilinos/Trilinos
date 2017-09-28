
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  directoryTest_Kokkos
  SOURCES directoryTest_Kokkos.cpp
  NUM_MPI_PROCS 4
  COMM mpi
)

# In progress to implement optimization such as functors.
# This was not brought to a working state yet.
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  directoryTest_KokkosOptimized
#  SOURCES directoryTest_KokkosOptimized.cpp
#  NUM_MPI_PROCS 4
#  COMM mpi
#)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  directoryTest_KokkosSimple
  SOURCES directoryTest_KokkosSimple.cpp
  NUM_MPI_PROCS 4
  COMM mpi
)

# Temporary tests for development and comparison
# with above Kokkos tests. To be deleted.

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


# These test are labeled PERFORMANCE and being
# used for development purposes.
# Perhaps all to be deleted eventually

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  directoryTest_PerformanceOriginal
  CATEGORIES PERFORMANCE
  SOURCES directoryTest_PerformanceOriginal.cpp
  NUM_MPI_PROCS 4
  COMM mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  directoryTest_PerformanceRelic
  CATEGORIES PERFORMANCE
  SOURCES directoryTest_PerformanceRelic.cpp
  NUM_MPI_PROCS 4
  COMM mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  directoryTest_PerformanceTpetra
  CATEGORIES PERFORMANCE
  SOURCES directoryTest_PerformanceTpetra.cpp
  NUM_MPI_PROCS 4
  COMM mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  directoryTest_PerformanceKokkos
  CATEGORIES PERFORMANCE
  SOURCES directoryTest_PerformanceKokkos.cpp
  NUM_MPI_PROCS 4
  COMM mpi
)

# In progress to implement optimization such as functors.
# This was not brought to a working state yet.
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  directoryTest_PerformanceKokkosOptimized
#  CATEGORIES PERFORMANCE
#  SOURCES directoryTest_PerformanceKokkosOptimized.cpp
#  NUM_MPI_PROCS 4
#  COMM mpi
#)
