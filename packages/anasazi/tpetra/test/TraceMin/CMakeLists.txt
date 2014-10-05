
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

# FIXME (mfh 25 Aug 2014) The "smallest" tests in this directory are
# running forever without converging.  Not sure why that is.  I'll
# work with Alicia to fix them, but for now I'll disable them.  Note
# that the TraceMin tests that look for the _largest_ eigenvalues pass
# just fine, so it must be a problem with shift and invert, or with
# the Belos solver in the tests.
#
# FIXME (mfh 15 Sep 2014) All the tests in this file now pass, except
# for the block-diagonal preconditioned MINRES test below.  That one
# runs for quite a long time.  I will leave it disabled for now.

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   Tpetra_TraceMin_smallest_bdprec_test
#   SOURCES cxx_main_standard_noprec.cpp
#   ARGS --saddleSolType="Block Diagonal Preconditioned Minres"
#   COMM serial mpi
#   )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_TraceMin_largest_standard_test
  SOURCES cxx_main_standard_noprec.cpp 
  ARGS --which=LM
  COMM serial mpi
  )
