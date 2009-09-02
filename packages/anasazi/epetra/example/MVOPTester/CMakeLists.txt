INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE(
  Epetra_MVOPTester_example
  SOURCES MVOPTesterEx.cpp
  COMM serial mpi
)
