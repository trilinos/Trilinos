
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  BLAS_test  
  SOURCES cxx_main.cpp
  ARGS -v
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

#PACKAGE_ADD_EXECUTABLE_AND_TEST(
#  BLAS_tmpl_test  
#  SOURCES cxx_tmpl_main.cpp
#  ARGS -v
#  NUM_MPI_PROCS 1
#  COMM serial mpi
#  STANDARD_PASS_OUTPUT
#  )
