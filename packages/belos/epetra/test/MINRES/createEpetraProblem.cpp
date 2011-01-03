//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#include "createEpetraProblem.hpp"
#include "Teuchos_Workspace.hpp"
#include "Trilinos_Util.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

int Belos::createEpetraProblem(
			       std::string                      &filename
			       ,RCP<Epetra_Map>         *rowMap
			       ,RCP<Epetra_CrsMatrix>   *A
			       ,RCP<Epetra_MultiVector> *B
			       ,RCP<Epetra_MultiVector> *X
			       ,int                             *MyPID_out
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
  epetraComm = rcp(new Epetra_MpiComm( MPI_COMM_WORLD ) );	
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
  *A = rcp(new Epetra_CrsMatrix(Copy, *epetraMap, NumNz));
  Teuchos::set_extra_data( epetraMap, "Operator::Map", Teuchos::ptr(A) );
  //
  // Add rows one-at-a-time
  //
  int NumEntries;
  for (i=0; i<NumMyElements; i++) {
    row_vals = val + bindx[i];
    col_inds = bindx + bindx[i];
    NumEntries = bindx[i+1] - bindx[i];
    int info = (*A)->InsertGlobalValues(update[i], NumEntries, row_vals, col_inds);
    assert( info == 0 );
    info = (*A)->InsertGlobalValues(update[i], 1, val+i, update+i);
    assert( info == 0 );
  }
  //
  // Finish up
  //
  int info = (*A)->FillComplete();
  assert( info == 0 );
  info = (*A)->OptimizeStorage();
  assert( info == 0 );
  (*A)->SetTracebackMode(1); // Shutdown Epetra Warning tracebacks
  //
  // Construct the right-hand side and solution multivectors.
  //
  if(B) {
    *B = rcp(new Epetra_MultiVector(::Copy, *epetraMap, b, NumMyElements, 1 ));
    Teuchos::set_extra_data( epetraMap, "B::Map", Teuchos::ptr(B) );
  }
  if(X) {
    *X = rcp(new Epetra_MultiVector(*epetraMap, 1 ));
    Teuchos::set_extra_data( epetraMap, "X::Map", Teuchos::ptr(X) );
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
