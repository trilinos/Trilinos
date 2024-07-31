// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// CrsGraph_BTF Test routine

#include <EpetraExt_ConfigDefs.h>
#include "EpetraExt_Version.h"

#include "Epetra_Time.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif

#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_LinearProblem.h"
#include "EpetraExt_AmesosBTF_CrsGraph.h"
#include "EpetraExt_LPTrans_From_GraphTrans.h"
#include "EpetraExt_Reindex_LinearProblem.h"

int main(int argc, char *argv[]) {

  int returnierr=0;

  using std::cout;
  using std::endl;
  using std::flush;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (size > 1) {
    cout << "This example cannot be run on more than one processor!" << endl;
    MPI_Finalize();
    returnierr = -1;
    return returnierr;
  }

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

  if (verbose) {
    cout << EpetraExt::EpetraExt_Version() << endl << endl;
    cout << Comm << endl << flush;
  }

  Comm.Barrier();

  int NumMyElements = 3;
 
  Epetra_Map Map( NumMyElements, 0, Comm );
  
  Epetra_CrsGraph Graph( Copy, Map, 1 );

  int index[2];
  index[0] = 2;
  Graph.InsertGlobalIndices( 0, 1, &index[0] );
  index[0] = 0;
  index[1] = 2;
  Graph.InsertGlobalIndices( 1, 2, &index[0] );
  index[0] = 1;
  Graph.InsertGlobalIndices( 2, 1, &index[0] );

  Graph.FillComplete();
  if (verbose) {
    cout << "***************** PERFORMING BTF TRANSFORM ON CRS_GRAPH *****************" <<endl<<endl;
    cout << "CrsGraph *before* BTF transform: " << endl << endl;
    cout << Graph << endl;
  }

  EpetraExt::AmesosBTF_CrsGraph BTFTrans( true, verbose );
  Epetra_CrsGraph & NewBTFGraph = BTFTrans( Graph );

  if (verbose) {
    cout << "CrsGraph *after* BTF transform: " << endl << endl;
    cout << NewBTFGraph << endl;
  }
	
  // Use BTF graph transformation to solve linear system.
  // Create an Epetra::CrsMatrix
  Epetra_CrsMatrix Matrix( Copy, Graph );
  double value[2];
  index[0] = 2; value[0] = 3.0;
  Matrix.ReplaceMyValues( 0, 1, &value[0], &index[0] );
  index[0] = 0; index[1] = 2;
  value[0] = 2.0; value[1] = 2.5;
  Matrix.ReplaceMyValues( 1, 2, &value[0], &index[0] );
  index[0] = 1; value[0] = 1.0;
  Matrix.ReplaceMyValues( 2, 1, &value[0], &index[0] );
  Matrix.FillComplete();

  // Create the solution and right-hand side vectors.
  Epetra_MultiVector RHS( Map, 1 ), LHS( Map, 1 );
  LHS.PutScalar( 0.0 );
  RHS.ReplaceMyValue( 0, 0, 3.0 );
  RHS.ReplaceMyValue( 1, 0, 4.5 );
  RHS.ReplaceMyValue( 2, 0, 1.0 );
  Epetra_LinearProblem problem( &Matrix, &LHS, &RHS );

  if (verbose) {
    cout << "*************** PERFORMING BTF TRANSFORM ON LINEAR_PROBLEM **************" <<endl<<endl;
    cout << "CrsMatrix *before* BTF transform: " << endl << endl;
    cout << Matrix << endl;
    cout << "MultiVector RHS *before* BTF transform: " << endl << endl;
    RHS.Print( cout );
  }

  // Create the linear problem transform.
  EpetraExt::LinearProblem_GraphTrans * LPTrans =
        new EpetraExt::LinearProblem_GraphTrans(
        *(dynamic_cast<EpetraExt::StructuralSameTypeTransform<Epetra_CrsGraph>*>(&BTFTrans)) );
  Epetra_LinearProblem* tProblem = &((*LPTrans)( problem ));
  LPTrans->fwd();

  if (verbose) {
    cout << "CrsMatrix *after* BTF transform: " << endl << endl;
    dynamic_cast<Epetra_CrsMatrix*>(tProblem->GetMatrix())->Print( cout );
    cout << "MultiVector RHS *after* BTF transform: " << endl << endl;
    tProblem->GetRHS()->Print( cout );
  }

  if (verbose) {
    cout << endl << "*************** PERFORMING REINDEXING ON LINEAR_PROBLEM **************" <<endl<<endl;
  }
  EpetraExt::ViewTransform<Epetra_LinearProblem> * ReIdx_LPTrans =
        new EpetraExt::LinearProblem_Reindex( 0 );

  Epetra_LinearProblem* tProblem2 = &((*ReIdx_LPTrans)( *tProblem ));
  ReIdx_LPTrans->fwd();

  if (verbose) {
    cout << endl << "CrsMatrix *after* BTF transform *and* reindexing: " << endl << endl;
    dynamic_cast<Epetra_CrsMatrix*>(tProblem2->GetMatrix())->Print( cout );
    cout << endl <<"Column Map *before* reindexing: " << endl << endl;
    cout << dynamic_cast<Epetra_CrsMatrix*>(tProblem->GetMatrix())->ColMap() << endl;
    cout << "Column Map *after* reindexing: " << endl << endl;
    cout << dynamic_cast<Epetra_CrsMatrix*>(tProblem2->GetMatrix())->ColMap() << endl;
  }
  
#ifdef EPETRA_MPI
  MPI_Finalize();
#endif
  
  return returnierr;
}
  
