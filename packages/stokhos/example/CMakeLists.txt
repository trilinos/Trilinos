INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  hermite_example
  SOURCES hermite_example.cpp
  ARGS -v
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  product_basis_example
  SOURCES product_basis_example.cpp
  ARGS -v
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE(
  dakota_pce
  SOURCES dakota_pce.cpp
  COMM serial mpi
  )
