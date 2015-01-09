//@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER

#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Ifpack_Reordering.h"
#include "Ifpack_RCMReordering.h"
#include "Ifpack_ReorderFilter.h"
#include "Ifpack_Utils.h"
#include "Teuchos_RefCountPtr.hpp"

//==============================================================================
int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // only one process
  if (Comm.NumProc() != 1) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    if (Comm.MyPID() == 0)
      cout << "Please run this test with one process only" << endl;
    // return success not to break the tests
    exit(EXIT_SUCCESS);
  }

  // ======================================================== //
  // now create the famous "upper arrow" matrix, which        //
  // should be reordered as a "lower arrow". Sparsity pattern //
  // will be printed on screen.                               //
  // ======================================================== //
  
  int NumPoints = 16;
  
#if !defined(EPETRA_NO_32BIT_GLOBAL_INDICES) || !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
  Epetra_Map Map(-1,NumPoints,0,Comm);
#else
  Epetra_Map Map;
#endif

  std::vector<int> Indices(NumPoints);
  std::vector<double> Values(NumPoints);

  Teuchos::RefCountPtr<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(Copy,Map,0) );
  for (int i = 0 ; i < NumPoints ; ++i) {
    
    int NumEntries;
    if (i == 0) {
      NumEntries = NumPoints;
      for (int j = 0 ; j < NumPoints ; ++j) {
	Indices[j] = j;
	Values[j] = 1.0;
      }
    }
    else {
      NumEntries = 2;
      Indices[0] = 0;
      Indices[1] = i;
      Values[0] = 1.0;
      Values[1] = 1.0;
    }

#if !defined(EPETRA_NO_32BIT_GLOBAL_INDICES) || !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
    A->InsertGlobalValues(i, NumEntries, &Values[0], &Indices[0]);
#endif
  }

  A->FillComplete();

  // print the sparsity to file, postscript format
  ////Ifpack_PrintSparsity(A,"OrigA.ps");

  // create the reordering...
  Teuchos::RefCountPtr<Ifpack_RCMReordering> Reorder = Teuchos::rcp( new Ifpack_RCMReordering() );
  // and compute is on A
  IFPACK_CHK_ERR(Reorder->Compute(*A));

  // cout information
  cout << *Reorder;

  // create a reordered matrix
  Ifpack_ReorderFilter ReordA(A, Reorder);

  // print the sparsity to file, postscript format
  ////Ifpack_PrintSparsity(ReordA,"ReordA.ps");

#ifdef HAVE_MPI
  MPI_Finalize(); 
#endif
  return(EXIT_SUCCESS);

}
