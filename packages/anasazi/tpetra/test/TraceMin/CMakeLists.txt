
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_TraceMin_smallest_proj_test
  SOURCES cxx_main_standard_noprec.cpp 
  ARGS 
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_TraceMin_smallest_schur_test
  SOURCES cxx_main_standard_noprec.cpp
  ARGS --saddleSolType="Schur Complement"
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_TraceMin_smallest_bdprec_test
  SOURCES cxx_main_standard_noprec.cpp
  ARGS --saddleSolType="Block Diagonal Preconditioned Minres"
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_TraceMin_largest_standard_test
  SOURCES cxx_main_standard_noprec.cpp 
  ARGS --which=LM
  COMM serial mpi
  )
