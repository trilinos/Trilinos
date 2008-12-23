INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_epetra_adapters
  SOURCES test_epetra_adapters.cpp
  ARGS --verbose
  COMM serial mpi
  )
