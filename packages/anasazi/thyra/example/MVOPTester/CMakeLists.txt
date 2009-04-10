INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE(
  Thyra_MVOPTester_example
  SOURCES MVOPTesterEx.cpp
  COMM serial mpi
)
