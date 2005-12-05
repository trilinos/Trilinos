#include <iostream>
#include <mpi.h>

using namespace std;


int
main(int argc, char *argv[])
{
  MPI::Init(argc, argv);
  
  cout << "Hello World! MPI works"  << endl ; 
  
  MPI::Finalize();
  return 0;
}
