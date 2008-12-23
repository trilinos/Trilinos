INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Map
  SOURCES Map.cpp
  COMM serial mpi
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix
  SOURCES CrsMatrix.cpp
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  VbrMatrix
  SOURCES VbrMatrix.cpp
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  LinearProblem
  SOURCES LinearProblem.cpp
  COMM serial mpi
  )

