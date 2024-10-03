// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "BelosEpetraUtils.h"
#include "BelosEpetraAdapter.hpp"
#include "Teuchos_Workspace.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "BelosGlobalComm.hpp"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

#ifdef HAVE_BELOS_TRIUTILS
#include "Trilinos_Util.h"
#endif

// RWH TODO FIXME: Trilinos_Util_distrib_msr_matrix available only if 32 bit GIDs available.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES

namespace Belos {

namespace Util {

#ifdef HAVE_BELOS_TRIUTILS

int createEpetraProblem( std::string              &filename,
                         RCP<Epetra_Map>          *rowMap,
                         RCP<Epetra_CrsMatrix>    *A,
                         RCP<Epetra_MultiVector>  *B,
                         RCP<Epetra_MultiVector>  *X,
                         int                      *MyPID
                       )
{
  int one = 1;
  return Belos::Util::createEpetraProblem( filename, rowMap, A, B, X, MyPID, one );
}

int createEpetraProblem( std::string             &filename,
			 RCP<Epetra_Map>         *rowMap,
			 RCP<Epetra_CrsMatrix>   *A,
			 RCP<Epetra_MultiVector> *B,
			 RCP<Epetra_MultiVector> *X,
			 int                     *MyPID_out,
                         int                     &numRHS
		       )
{
  //
  int &MyPID = *MyPID_out;
  //
  int i;
  int n_nonzeros, N_update;
  int *bindx=0, *update=0, *col_inds=0;
  double *val=0, *row_vals=0;
  double *xguess=0, *b=0, *xexact=0;

  RCP<Epetra_Comm> epetraComm;
#ifdef EPETRA_MPI	
  epetraComm = rcp(new Epetra_MpiComm( Belos::get_global_comm() ) );
#else	
  epetraComm = rcp(new Epetra_SerialComm());
#endif
	
  MyPID = epetraComm->MyPID();
  //
  // **********************************************************************
  // ******************Set up the problem to be solved*********************
  // **********************************************************************
  //
  int NumGlobalElements;  // total # of rows in matrix
  //
  // *****Read in matrix from HB file******
  //
  Trilinos_Util_read_hb(const_cast<char *>(filename.c_str()), MyPID, &NumGlobalElements, &n_nonzeros,
			&val, &bindx, &xguess, &b, &xexact);
  // 
  // *****Distribute data among processors*****
  //
  Trilinos_Util_distrib_msr_matrix(*epetraComm, &NumGlobalElements, &n_nonzeros, &N_update, 
				   &update, &val, &bindx, &xguess, &b, &xexact);
  //
  // *****Construct the matrix*****
  //
  int NumMyElements = N_update; // # local rows of matrix on processor
  //
  // Create an integer std::vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
  // on this processor
  //
  int * NumNz = new int[NumMyElements];
  for (i=0; i<NumMyElements; i++) {
    NumNz[i] = bindx[i+1] - bindx[i] + 1;
  }
  //
  RCP<Epetra_Map> epetraMap = rcp(new Epetra_Map(NumGlobalElements, NumMyElements, update, 0, *epetraComm));
  Teuchos::set_extra_data( epetraComm, "Map::Comm", Teuchos::inOutArg(epetraMap) );
  if(rowMap) *rowMap = epetraMap;
  //
  // Create a Epetra_Matrix
  //
  *A = rcp(new Epetra_CrsMatrix(BELOSEPETRACOPY, *epetraMap, NumNz));
  Teuchos::set_extra_data( epetraMap, "Operator::Map", Teuchos::ptr(A) );
  //
  // Add rows one-at-a-time
  //
  int NumEntries;
  for (i=0; i<NumMyElements; i++) {
    row_vals = val + bindx[i];
    col_inds = bindx + bindx[i];
    NumEntries = bindx[i+1] - bindx[i];
#ifndef NDEBUG
    int info =
#endif
    (*A)->InsertGlobalValues(update[i], NumEntries, row_vals, col_inds);
    assert( info == 0 );
#ifndef NDEBUG
    info =
#endif
    (*A)->InsertGlobalValues(update[i], 1, val+i, update+i);
    assert( info == 0 );
  }
  //
  // Finish up
  //
#ifndef NDEBUG
  int info =
#endif
  (*A)->FillComplete();
  assert( info == 0 );
#ifndef NDEBUG
  info =
#endif
  (*A)->OptimizeStorage();
  assert( info == 0 );
  (*A)->SetTracebackMode(1); // Shutdown Epetra Warning tracebacks
  //
  // Construct the right-hand side and solution multivectors.
  //
  if (B) {
    if (b != NULL) {
      *B = rcp(new Epetra_MultiVector(BELOSEPETRACOPY, *epetraMap, b, NumMyElements, 1 ));
      Teuchos::set_extra_data( epetraMap, "B::Map", Teuchos::ptr(B) );
      numRHS = 1;
    }
    else {
      *B = rcp (new Epetra_MultiVector (*epetraMap, numRHS));
      Teuchos::set_extra_data( epetraMap, "B::Map", Teuchos::ptr(B) );
      (*B)->Random ();
    }
  }
  if (X) {
    *X = rcp (new Epetra_MultiVector (*epetraMap, numRHS));
    Teuchos::set_extra_data( epetraMap, "X::Map", Teuchos::ptr(X) );
    (*X)->PutScalar (0.0);
  } 
  //
  // Create workspace
  //
  Teuchos::set_default_workspace_store(
    Teuchos::rcp(new Teuchos::WorkspaceStoreInitializeable(static_cast<size_t>(2e+6)))
    );
  //
  // Free up memory
  //
  delete [] NumNz;
  free(update);
  free(val);
  free(bindx);
  if (xexact) free(xexact);
  if (xguess) free(xguess);
  if (b) free(b);

  return (0);
}

#endif

int rebalanceEpetraProblem( RCP<Epetra_Map>         &Map,
                            RCP<Epetra_CrsMatrix>   &A,
                            RCP<Epetra_MultiVector> &B,
                            RCP<Epetra_MultiVector> &X,
                            Epetra_Comm             &Comm
                          )
{
  // Rebalance linear system across multiple processors.
  if ( Comm.NumProc() > 1 ) {
    RCP<Epetra_Map> newMap = rcp( new Epetra_Map( Map->NumGlobalElements(), Map->IndexBase(), Comm ) );
    RCP<Epetra_Import> newImport = rcp( new Epetra_Import( *newMap, *Map ) );

    // Create rebalanced versions of the linear system.
    RCP<Epetra_CrsMatrix> newA = rcp( new Epetra_CrsMatrix( BELOSEPETRACOPY, *newMap, 0 ) );
    newA->Import( *A, *newImport, Insert );
    newA->FillComplete();
    RCP<Epetra_MultiVector> newB = rcp( new Epetra_MultiVector( *newMap, B->NumVectors() ) );
    newB->Import( *B, *newImport, Insert );
    RCP<Epetra_MultiVector> newX = rcp( new Epetra_MultiVector( *newMap, X->NumVectors() ) );
    newX->Import( *X, *newImport, Insert );

    // Set the pointers to the new rebalance linear system.
    A = newA;
    B = newB;
    X = newX;
    Map = newMap;
  }

  return (0);
}

} // namespace Util

namespace Test {

    MPISession::MPISession (Teuchos::Ptr<int> argc, Teuchos::Ptr<char**> argv)
    {
#ifdef EPETRA_MPI
      MPI_Init (argc.getRawPtr(), argv.getRawPtr());
#endif // EPETRA_MPI
    }

    MPISession::~MPISession ()
    {
#ifdef EPETRA_MPI
      MPI_Finalize ();
#endif // EPETRA_MPI
    }

    Teuchos::RCP<const Epetra_Comm>
    MPISession::getComm ()
    {
      using Teuchos::rcp;

      if (comm_.is_null()) {
#ifdef EPETRA_MPI
        comm_ = rcp (new Epetra_MpiComm (Belos::get_global_comm()));
#else
        comm_ = rcp (new Epetra_SerialComm);
#endif // EPETRA_MPI
      }
      return comm_;
    }

} // namespace Test

} // namespace Belos

#endif // EPETRA_NO_32BIT_GLOBAL_INDICES
