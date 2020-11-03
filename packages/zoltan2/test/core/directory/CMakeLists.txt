
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  directoryTest_Kokkos
  SOURCES directoryTest_Kokkos.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  directoryTest_KokkosSimple
  SOURCES directoryTest_KokkosSimple.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
)

# The purpose of this test is to duplicate the
# result of findUniqueGids but using the new
# Zoltan2Directory class.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  directoryTest_findUniqueGids.cpp
  SOURCES directoryTest_findUniqueGids.cpp
  COMM serial mpi
  FAIL_REGULAR_EXPRESSION "FAIL"
)
