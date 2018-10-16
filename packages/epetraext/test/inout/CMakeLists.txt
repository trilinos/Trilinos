

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  inout_test
  SOURCES
    cxx_main.cpp
    Poisson2dOperator.cpp
  ARGS -v
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  inout_hdf5_test
  SOURCES
    EpetraExt_HDF5_UnitTest.cpp
  COMM serial mpi
  )
