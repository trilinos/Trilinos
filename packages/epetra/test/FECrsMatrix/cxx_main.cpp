// Epetra_FECrsMatrix Test routine

#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_FECrsMatrix.h"
#include "ExecuteTestProblems.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "../epetra_test_err.h"

int main(int argc, char *argv[]) {

  int ierr = 0;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;
  Epetra_SerialComm Comm;

#endif

//  Comm.SetTracebackMode(0); // This should shut down any error tracing
  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

#ifdef EPETRA_MPI
  int localverbose = verbose ? 1 : 0;
  int globalverbose=0;
  MPI_Allreduce(&localverbose, &globalverbose, 1, MPI_INT, MPI_SUM,
		MPI_COMM_WORLD);
  verbose = (globalverbose>0);
#endif

  //  char tmp;
  //  if (rank==0) cout << "Press any key to continue..."<< endl;
  //  if (rank==0) cin >> tmp;
  //  Comm.Barrier();

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc(); 
  if (verbose) cout << Comm <<endl;

  // Redefine verbose to only print on PE 0
  //if (verbose && rank!=0) verbose = false;

  int NumVectors = 1;
  int NumMyElements = 4;
  int NumGlobalElements = NumMyElements*NumProc;
  int IndexBase = 0;
  
  // Test Petra-defined uniform linear distribution constructor

  if (verbose&&rank==0) {
    cout << "\n*********************************************************"<<endl;
    cout << "Checking Epetra_Map(NumGlobalElements, NumMyElements, IndexBase, Comm)" << endl;
    cout << "*********************************************************" << endl;
  }

  Epetra_Map Map(NumGlobalElements, NumMyElements, IndexBase, Comm);

  EPETRA_TEST_ERR( Drumm1(Map, verbose),ierr);

  EPETRA_TEST_ERR( Drumm2(Map, verbose),ierr);

  EPETRA_TEST_ERR( Drumm3(Map, verbose),ierr);

  EPETRA_TEST_ERR(MultiVectorTests(Map, NumVectors, verbose),ierr);

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return ierr;
}

