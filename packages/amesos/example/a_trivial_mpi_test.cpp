#include <iostream>
#include <mpi.h>

using namespace std;


int
main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);  // MPI_Init(argc,argv) won't compile on IRIX64 - see bug #1885
  
  cout << "Hello World! MPI works"  << endl ; 
  
  MPI_Finalize();
  return 0;
}
