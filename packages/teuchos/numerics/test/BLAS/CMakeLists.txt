
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BLAS_test  
  SOURCES cxx_main.cpp
  ARGS -v
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BLAS_tmpl_test  
  SOURCES cxx_tmpl_main.cpp
  ARGS -v
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BLAS_tmpl_comp_test  
  SOURCES cxx_tmpl_main_comp.cpp
  ARGS -v
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BLAS_ROTG_test
  SOURCES default_blas_rot.cpp
  ARGS "--verbose"
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
