// Tpetra Platform tester
// Modified: 21-Jan-2003

#define SCALARTYPE float
#define ORDINALTYPE int

#ifdef TPETRA_MPI
#include "Tpetra_MpiPlatform.hpp" 
#include <mpi.h>
#else
#include "Tpetra_SerialPlatform.hpp"
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
  Tpetra::MpiPlatform<SCALARTYPE, ORDINALTYPE> platform( MPI_COMM_WORLD );
#else
  int size = 1; // Serial case (not using MPI)
  int rank = 0;
  if(verbose) cout << "Creating SerialPlatform object...";
  Tpetra::SerialPlatform<SCALARTYPE, ORDINALTYPE> platform;
  //if(debug) cout << platform.label() << endl;
	if(verbose) cout << "Successful." << endl;
#endif

	if(verbose) cout << "Creating Comm objects...";
	Tpetra::Comm<SCALARTYPE, ORDINALTYPE>* comm1 = platform.createScalarComm();
	Tpetra::Comm<ORDINALTYPE, ORDINALTYPE>* comm2 = platform.createOrdinalComm();
	delete comm1;
	delete comm2;
	if(verbose) cout << "Successful." << endl;

	if(verbose) cout << "Creating Distributor objects...";
	Tpetra::Distributor<SCALARTYPE, ORDINALTYPE>* distributor1 = platform.createScalarDistributor();
	Tpetra::Distributor<ORDINALTYPE, ORDINALTYPE>* distributor2 = platform.createOrdinalDistributor();
	delete distributor1;
	delete distributor2;
	if(verbose) cout << "Successful." << endl;
  
	cout << "Platform test successful." << endl;
  
  return(0);
}
