// Tpetra Platform tester
// Modified: 14-Oct-2002

#define PACKETTYPE int
#define ORDINALTYPE int

#ifdef TPETRA_MPI
#include "Tpetra_MpiPlatform.h" 
#include <mpi.h>
#else
#include "Tpetra_SerialPlatform.h"
#endif

int main(int argc, char* argv[]) {
	bool verbose = false;
	if (argc>1 && argv[1][0]=='-' && argv[1][1]=='v') 
		verbose = true;

#ifdef TPETRA_MPI 
  // Initialize MPI
  MPI_Init(&argc,&argv);complex<float>
  int size, rank; // Number of MPI processes, My process ID
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Tpetra::MpiPlatform<PACKETTYPE, ORDINALTYPE> platform( MPI_COMM_WORLD );
#else
  int size = 1; // Serial case (not using MPI)
  int rank = 0;
  Tpetra::SerialPlatform<PACKETTYPE, ORDINALTYPE> platform;
#endif

  if(verbose) cout << "\nPlatform object created." << endl;
  if(verbose) cout << platform.label() << endl;
  
  if(verbose) cout << "Testing getImageID and getNumImages...";
  assert(platform.getMyImageID() == rank);
  assert(platform.getNumImages() == size);
	if(verbose) cout << "Successful." << endl;

	if(verbose) cout << "Creating Comm object...";
	Tpetra::Comm<PACKETTYPE, ORDINALTYPE>* comm = platform.createComm();
	delete comm;
	if(verbose) cout << "Successful." << endl;
  
	cout << "Platform test successful." << endl;
  
  return(0);
}
