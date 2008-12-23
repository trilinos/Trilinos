
INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MultiLevelPreconditioner_NonSym
  SOURCES MultiLevelPreconditioner_NonSym.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )
