// CrsGraph_BTF Test routine

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "ET_CrsGraph_BTF.h"
#include "../epetra_test_err.h"

int main(int argc, char *argv[]) {

  int i, ierr=0, returnierr=0;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;


#ifdef EPETRA_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  if (!verbose) Comm.SetTracebackMode(0); // This should shut down any error traceback reporting

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  if (verbose) cout << Comm << endl << flush;

  Comm.Barrier();
  bool verbose1 = verbose;

  if (verbose) verbose = (MyPID==0);

  int NumMyElements = 3;
  int NumGlobalElements = NumMyElements;
  int IndexBase = 0;
  
  Epetra_Map Map( NumMyElements, 0, Comm );
  cout << Map << endl;
  
  Epetra_CrsGraph Graph( Copy, Map, 1 );

  int index = 2;
  Graph.InsertGlobalIndices( 0, 1, &index );
  index = 0;
  Graph.InsertGlobalIndices( 1, 1, &index );
  index = 1;
  Graph.InsertGlobalIndices( 2, 1, &index );

  Graph.TransformToLocal();
  cout << Graph << endl;

  Epetra_Transform::SameTypeTransform<Epetra_CrsGraph> * BTF_Trans
		= new Epetra_Transform::CrsGraph_BTF();

  std::auto_ptr<Epetra_CrsGraph> NewGraph( (*BTF_Trans)( Graph ) );

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return returnierr;
}

