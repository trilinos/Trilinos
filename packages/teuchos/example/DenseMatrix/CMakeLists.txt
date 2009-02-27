
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  DenseMatrix_example
  COMM serial mpi  
  NUM_MPI_PROCS 1
  SOURCES cxx_main.cpp
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  SymDenseMatrix_example
  COMM serial mpi  
  NUM_MPI_PROCS 1
  SOURCES cxx_main_sym.cpp
  )
