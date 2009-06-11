
INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Test_Epetra_CrsMatrix
  SOURCES cxx_main.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  )
