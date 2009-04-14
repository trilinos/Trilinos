INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Map
  SOURCES Map.cpp
  NUM_MPI_PROCS 1
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix
  SOURCES CrsMatrix.cpp
  NUM_MPI_PROCS 1
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  VbrMatrix
  SOURCES VbrMatrix.cpp
  NUM_MPI_PROCS 1
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  LinearProblem
  SOURCES LinearProblem.cpp
  NUM_MPI_PROCS 1
  COMM serial mpi
  )

