//  See "copypright.h" for copyright information
#include "copyright.h"
//
// hello_test
//
//  usage: 
//     hello_test
//
//  output:  
//     prints a summary line and one line "Hello" for each process to standard out
//     If --enable-newp_swahili is set on the configure line:
//        prints a summary line and one line "Jambo" for each process to standard out
//
#endif
#include "Newp_Hello.h"
#ifdef HAVE_NEWP_SWAHILI
#include "Newp_Jambo.h"
#endif

int main(int argc, char **argv)
{
  //
  //  If --enable-mpi, an MPI communicator is used, otherwise a serial
  //  stub communicator is used.  
  //
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  //
  //  Print out a summary line followed by a "Hello" line from each process
  //
  Newp_Hello Hello( Comm ) ; 
  Hello.Print( cout );


  //
  //  If --enable-newp_swahili is set, HAVE_NEWP_SWAHILI is set in 
  //    New_Package_config.h which is included by Newp_Hello.h and hence:
  //      Print out a summary line followed by a "Jambo" line from each process
  //
#ifdef HAVE_NEWP_SWAHILI
  Newp_Jambo Jambo( Comm ) ; 
  Jambo.Print( cout );
#endif

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif
  return 0;
}

  

