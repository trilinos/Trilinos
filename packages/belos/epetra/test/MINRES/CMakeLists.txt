
INCLUDE(PackageAddExecutableAndTest)

# FIXME (mfh 06 Dec 2010) Nico's MINRES test executable is broken when
# running on >1 MPI process.  Our temporary solution is to limit
# execution of the MPI test to 1 MPI process.
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Epetra_MINRES
  SOURCES test_minres.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS
    "--debug --verbose --tol=1e-4"
  STANDARD_PASS_OUTPUT 
  )
