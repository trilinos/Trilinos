// Transform_Composite Test routine

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif

#include "Epetra_SerialComm.h"
#include "Epetra_Time.h"

#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"

#include "Epetra_Transform_Composite.h"
#include "EDT_LinearProblem_GraphTrans.h"
#include "EDT_CrsGraph_SymmRCM.h"
#include "EDT_CrsGraph_Overlap.h"
#include "EDT_CrsGraph_BTF.h"

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

  Epetra_Map ViewMap( NumMyElements-2, 0, Comm );
  cout << ViewMap << endl;
  
  Epetra_CrsGraph Graph( Copy, Map, 1 );

  int index = 2;
  Graph.InsertGlobalIndices( 0, 1, &index );
  index = 0;
  Graph.InsertGlobalIndices( 1, 1, &index );
  index = 1;
  Graph.InsertGlobalIndices( 2, 1, &index );

  Graph.TransformToLocal();
  cout << Graph << endl;

  EpetraExt::Transform_Composite<Epetra_LinearProblem> CompTrans;

  EpetraExt::CrsGraph_BTF BTF_Trans;
  EpetraExt::SameTypeTransform<Epetra_LinearProblem> *
        BTF_LPTrans = new EpetraExt::LinearProblem_GraphTrans( BTF_Trans );
  CompTrans.addTransform( BTF_LPTrans );

  EpetraExt::CrsGraph_SymmRCM RCM_Trans;
  EpetraExt::SameTypeTransform<Epetra_LinearProblem> *
        RCM_LPTrans = new EpetraExt::LinearProblem_GraphTrans( RCM_Trans );
  CompTrans.addTransform( RCM_LPTrans );

#ifdef EPETRA_MPI
  EpetraExt::CrsGraph_Overlap Overlap_Trans(1);
  EpetraExt::SameTypeTransform<Epetra_LinearProblem> *
        Overlap_LPTrans = new EpetraExt::LinearProblem_GraphTrans( Overlap_Trans );
  CompTrans.addTransform( Overlap_LPTrans );
#endif

  Epetra_CrsMatrix Matrix( Copy, Graph );
  index = 2;
  double val = 2;
  Matrix.InsertGlobalValues( 0, 1, &val, &index );
  index = 0;
  val = 0;
  Matrix.InsertGlobalValues( 1, 1, &val, &index );
  index = 1;
  val = 1;
  Matrix.InsertGlobalValues( 2, 1, &val, &index);

  vector<double> valA(3);
  valA[0]=0; valA[1]=1; valA[2]=2;
  Epetra_BlockMap & MapRef = Map;
  Epetra_Vector LHS( Copy, MapRef, &valA[0] );
  Epetra_Vector RHS( Copy, MapRef, &valA[0] );

  Epetra_LinearProblem Prob( &Matrix, &LHS, &RHS );

  Epetra_LinearProblem & NewProb = CompTrans( Prob );

  CompTrans.fwd();
  CompTrans.rvs();

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return ierr;
}

