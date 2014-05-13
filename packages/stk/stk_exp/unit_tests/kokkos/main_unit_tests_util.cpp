#include <gtest/gtest.h>
#include <mpi.h>

int* STKUNIT_ARGC;
char** STKUNIT_ARGV;

int main( int argc, char **argv )
{
  ::testing::InitGoogleTest(&argc,argv);
  MPI_Init( & argc , & argv );
  STKUNIT_ARGC = &argc;
  STKUNIT_ARGV = argv;
  int err = RUN_ALL_TESTS();
  MPI_Finalize();
  return err;
}
