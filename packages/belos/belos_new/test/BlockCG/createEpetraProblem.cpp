// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "createEpetraProblem.hpp"
#include "Teuchos_Workspace.hpp"
#include "Trilinos_Util.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

int Belos::createEpetraProblem(
			 int                              argc
			 ,char                            *argv[]
			 ,RefCountPtr<Epetra_Map>         *rowMap
			 ,RefCountPtr<Epetra_CrsMatrix>   *A
			 ,RefCountPtr<Epetra_MultiVector> *B
			 ,RefCountPtr<Epetra_MultiVector> *X
			 ,int                             *MyPID_out
			 ,bool                            *verbose_out
			 )
{
  //
  int &MyPID = *MyPID_out;
  bool &verbose = *verbose_out;
  //
  int i;
  int n_nonzeros, N_update;
  int *bindx=0, *update=0, *col_inds=0;
  double *val=0, *row_vals=0;
  double *xguess=0, *b=0, *xexact=0;

  RefCountPtr<Epetra_Comm> epetraComm;
#ifdef EPETRA_MPI	
  epetraComm = rcp(new Epetra_MpiComm( MPI_COMM_WORLD ) );	
#else	
  epetraComm = rcp(new Epetra_SerialComm());
#endif
	
  MyPID = epetraComm->MyPID();
	
  verbose = 0;
  //
  if((argc < 2 || argc > 4)&& MyPID==0) {
    cerr << "Usage: " << argv[0]
         << " [ -v ] [ HB_filename ]" << endl
         << "where:" << endl
         << "-v                 - run test in verbose mode" << endl
         << "HB_filename        - filename and path of a Harwell-Boeing data set" << endl
         << endl;
    return(1);
  }
  //
  // Find verbosity flag
  //
  int file_arg = 1;
  for(i = 1; i < argc; i++)
    {
      if(argv[i][0] == '-' && argv[i][1] == 'v') {
	verbose = (MyPID == 0);
	if(i==1) file_arg = 2;
      }
    }
  //
  // **********************************************************************
  // ******************Set up the problem to be solved*********************
  // **********************************************************************
  //
  int NumGlobalElements;  // total # of rows in matrix
  //
  // *****Read in matrix from HB file******
  //
  Trilinos_Util_read_hb(argv[file_arg], MyPID, &NumGlobalElements, &n_nonzeros,
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
  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
  // on this processor
  //
  int * NumNz = new int[NumMyElements];
  for (i=0; i<NumMyElements; i++) {
    NumNz[i] = bindx[i+1] - bindx[i] + 1;
  }
  //
  RefCountPtr<Epetra_Map> epetraMap = rcp(new Epetra_Map(NumGlobalElements, NumMyElements, update, 0, *epetraComm));
  Teuchos::set_extra_data( epetraComm, "Map::Comm", &epetraMap );
  if(rowMap) *rowMap = epetraMap;
  //
  // Create a Epetra_Matrix
  //
  *A = rcp(new Epetra_CrsMatrix(Copy, *epetraMap, NumNz));
  Teuchos::set_extra_data( epetraMap, "Operator::Map", A );
  //
  // Add rows one-at-a-time
  //
  int NumEntries;
  for (i=0; i<NumMyElements; i++) {
    row_vals = val + bindx[i];
    col_inds = bindx + bindx[i];
    NumEntries = bindx[i+1] - bindx[i];
    assert((*A)->InsertGlobalValues(update[i], NumEntries, row_vals, col_inds)==0);
    assert((*A)->InsertGlobalValues(update[i], 1, val+i, update+i)==0);
  }
  //
  // Finish up
  //
  assert((*A)->TransformToLocal()==0);
  assert((*A)->OptimizeStorage()==0);
  (*A)->SetTracebackMode(1); // Shutdown Epetra Warning tracebacks
  //
  // Construct the right-hand side and solution multivectors.
  //
  if(B) {
    *B = rcp(new Epetra_MultiVector(::Copy, *epetraMap, b, NumMyElements, 1 ));
    Teuchos::set_extra_data( epetraMap, "B::Map", B );
  }
  if(X) {
    *X = rcp(new Epetra_MultiVector(*epetraMap, 1 ));
    Teuchos::set_extra_data( epetraMap, "X::Map", X );
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
