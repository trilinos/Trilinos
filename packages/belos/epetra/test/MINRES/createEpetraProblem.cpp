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
#include "Epetra_Comm.h"
#ifdef EPETRA_MPI
#  include "mpi.h"
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif // EPETRA_MPI
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_StandardCatchMacros.hpp"


namespace Belos {

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
	comm_ = rcp (new Epetra_MpiComm (MPI_COMM_WORLD));	
#else	
	comm_ = rcp (new Epetra_SerialComm);
#endif // EPETRA_MPI
      }
      return comm_;
    }
  } // namespace Test


void
createEpetraProblem (const Teuchos::RCP<const Epetra_Comm>& epetraComm,
		     const std::string& filename,
		     Teuchos::RCP<Epetra_Map>& rowMap,
		     Teuchos::RCP<Epetra_CrsMatrix>& A,
		     Teuchos::RCP<Epetra_MultiVector>& B,
		     Teuchos::RCP<Epetra_MultiVector>& X,
		     int& numRHS)
{
  using Teuchos::inOutArg;
  using Teuchos::ptr;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::set_extra_data;

  const int MyPID = epetraComm->MyPID();

  int n_nonzeros, N_update;
  int *bindx = NULL, *update = NULL, *col_inds = NULL;
  double *val = NULL, *row_vals = NULL;
  double *xguess = NULL, *b = NULL, *xexact = NULL;

  //
  // Set up the problem to be solved
  //
  int NumGlobalElements;  // total # of rows in matrix

  try {
    // Read in matrix from HB file
    Trilinos_Util_read_hb (const_cast<char *> (filename.c_str()), MyPID, 
			   &NumGlobalElements, &n_nonzeros, &val, &bindx, 
			   &xguess, &b, &xexact);

    // Distribute data among processors
    Trilinos_Util_distrib_msr_matrix (*epetraComm, &NumGlobalElements, 
				      &n_nonzeros, &N_update, &update, &val, 
				      &bindx, &xguess, &b, &xexact);
    //
    // Construct the matrix
    //
    int NumMyElements = N_update; // # local rows of matrix on processor
    //
    // Create an int array NumNz that is used to build the Petra Matrix.
    // NumNz[i] is the number of OFF-DIAGONAL terms for the i-th global
    // equation on this processor.
    //
    std::vector<int> NumNz (NumMyElements);
    for (int i = 0; i < NumMyElements; ++i) {
      NumNz[i] = bindx[i+1] - bindx[i] + 1;
    }
    rowMap = rcp (new Epetra_Map (NumGlobalElements, NumMyElements, 
				  update, 0, *epetraComm));
    //set_extra_data (epetraComm, "Map::Comm", inOutArg (rowMap));

    // Create an Epetra sparse matrix.
    if (NumMyElements == 0) {
      A = rcp (new Epetra_CrsMatrix (Epetra_DataAccess::Copy, *rowMap, static_cast<int*>(NULL)));
    } else {
      A = rcp (new Epetra_CrsMatrix (Epetra_DataAccess::Copy, *rowMap, &NumNz[0]));
    }
    //set_extra_data (rowMap, "Operator::Map", ptr (A));

    // Add rows to the sparse matrix one at a time.
    int NumEntries;
    for (int i = 0; i < NumMyElements; ++i) {
      row_vals = val + bindx[i];
      col_inds = bindx + bindx[i];
      NumEntries = bindx[i+1] - bindx[i];
      int info = A->InsertGlobalValues (update[i], NumEntries, row_vals, col_inds);
      TEUCHOS_TEST_FOR_EXCEPTION( info != 0, std::logic_error,
				  "Failed to insert global value into A." );
      info = A->InsertGlobalValues (update[i], 1, val+i, update+i);
      TEUCHOS_TEST_FOR_EXCEPTION( info != 0, std::logic_error,
				  "Failed to insert global value into A." );
    }
    // Finish initializing the sparse matrix.
    int info = A->FillComplete();
    TEUCHOS_TEST_FOR_EXCEPTION( info != 0, std::logic_error,
				"FillComplete() failed on the sparse matrix A.");
    info = A->OptimizeStorage();
    TEUCHOS_TEST_FOR_EXCEPTION( info != 0, std::logic_error,
				"OptimizeStorage() failed on the sparse matrix A.");
    A->SetTracebackMode (1); // Shut down Epetra warning tracebacks
    //
    // Construct the right-hand side and solution multivectors.
    //
    if (false && b != NULL) {
      B = rcp (new Epetra_MultiVector (Epetra_DataAccess::Copy, *rowMap, b, NumMyElements, 1));
      numRHS = 1;
    } else {
      B = rcp (new Epetra_MultiVector (*rowMap, numRHS));
      B->Random ();
    }
    X = rcp (new Epetra_MultiVector (*rowMap, numRHS));
    X->PutScalar (0.0);

    //set_extra_data (rowMap, "X::Map", Teuchos::ptr (X));
    //set_extra_data (rowMap, "B::Map", Teuchos::ptr (B));
  } 
  catch (std::exception& e) {
    // Free up memory before rethrowing the exception.
    if (update) free(update);
    if (val) free(val);
    if (bindx) free(bindx);
    if (xexact) free(xexact);
    if (xguess) free(xguess);
    if (b) free(b);
    throw e; // Rethrow the exception.
  } 
  catch (...) { // Epetra sometimes throws an int "object."
    // Free up memory before rethrowing the exception.
    if (update) free(update);
    if (val) free(val);
    if (bindx) free(bindx);
    if (xexact) free(xexact);
    if (xguess) free(xguess);
    if (b) free(b);
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
			       "Generating the Epetra problem to solve failed "
			       "with an unknown error.");
  }

  //if (false) {
    //
    // Create workspace
    //
    //Teuchos::set_default_workspace_store (rcp (new Teuchos::WorkspaceStoreInitializeable(static_cast<size_t>(2e+6))));
  //}

  //
  // Free up memory.
  //
  if (update) free(update);
  if (val) free(val);
  if (bindx) free(bindx);
  if (xexact) free(xexact);
  if (xguess) free(xguess);
  if (b) free(b);
}

} // namespace Belos
