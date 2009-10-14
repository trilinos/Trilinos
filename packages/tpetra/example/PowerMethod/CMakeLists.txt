INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  PowerMethodExample
  SOURCES Tpetra_Power_Method
  ARGS 1000
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)
