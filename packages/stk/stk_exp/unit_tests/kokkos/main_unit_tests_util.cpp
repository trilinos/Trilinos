#include <gtest/gtest.h>
#include <mpi.h>

int gl_argc;
char** gl_argv;

int main( int argc, char **argv )
{
  ::testing::InitGoogleTest(&argc,argv);
  MPI_Init( & argc , & argv );
  gl_argc = argc;
  gl_argv = argv;
  int err = RUN_ALL_TESTS();
  MPI_Finalize();
  return err;
}
