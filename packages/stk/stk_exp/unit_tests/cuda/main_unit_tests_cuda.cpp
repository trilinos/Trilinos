#include <gtest/gtest.h>
#include <stk_util/stk_config.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

int gl_argc;
char** gl_argv;

int main( int argc, char **argv )
{
  ::testing::InitGoogleTest(&argc,argv);
#ifdef HAVE_MPI
  MPI_Init( & argc , & argv );
#endif
  gl_argc = argc;
  gl_argv = argv;
  int err = RUN_ALL_TESTS();
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return err;
}
