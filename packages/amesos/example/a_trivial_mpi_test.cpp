#include <iostream>
#include <mpi.h>

int
main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);  // MPI_Init(argc,argv) won't compile on IRIX64 - see bug #1885
  
  std::cout << "Hello World! MPI works"  << std::endl ; 
  
  MPI_Finalize();
  return 0;
}
