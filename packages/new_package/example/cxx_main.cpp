//
// hello_test
//
//  usage: 
//     hello_test
//
//  output:  
//     stdout
//
//  exits with 0 if test completed (does not imply that the test passed)
//  exits with -1 if command line options or file permissions are wrong 
//

//  #include "New_Package_config.h" - I suspect that Epatra_config.h has everything we need 
//  and hence obviates the need for New_Package_config.h
//  Or maybe not - we might still need HAVE_NEWP_SWAHILI 
//  
#include "Newp_Hello.h"
#ifdef HAVE_NEWP_SWAHILI
#include "Newp_Jambo.h"
#endif

main(int argc, char **argv)
{
  int exit_value = 0 ; 

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  Newp_Hello Hello( Comm ) ; 
  Hello.Print( cout );
#ifdef HAVE_NEWP_SWAHILI
  Newp_Jambo Jambo( Comm ) ; 
  Jambo.Print( cout );
#endif

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif
  exit( exit_value ) ; 
}

  

