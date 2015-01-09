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
#include "Ifpack_DropFilter.h"
#include "Ifpack_SparsityFilter.h"
#include "Ifpack_SingletonFilter.h"
#include "Ifpack_Utils.h"
#include "Teuchos_RefCountPtr.hpp"

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  if (Comm.NumProc() != 1) {
    cerr << "This example must be run with one process only." << endl;
    // exit with success not to break the test harness
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_SUCCESS);
  }
  
  int NumPoints = 5;
#if !defined(EPETRA_NO_32BIT_GLOBAL_INDICES) || !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
  Epetra_Map Map(NumPoints,0,Comm);
#else
  Epetra_Map Map;
#endif

  Teuchos::RefCountPtr<Epetra_CrsMatrix> Matrix = Teuchos::rcp( new Epetra_CrsMatrix(Copy,Map,0) );

  std::vector<int> Indices(NumPoints);
  std::vector<double> Values(NumPoints);
  double Diag = 0.0;

  for (int i = 0 ; i < NumPoints ; ++i) {
    // add a diagonal
    Matrix->InsertGlobalValues(i,1,&Diag,&i);

    // add off-diagonals
    int NumEntries = 0;
    for (int j = i + 1 ; j < NumPoints ; ++j) {
      Indices[NumEntries] = j;
      Values[NumEntries] = 1.0 * (j - i);
      ++NumEntries;
    }
#if !defined(EPETRA_NO_32BIT_GLOBAL_INDICES) || !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
    Matrix->InsertGlobalValues(i,NumEntries,&Values[0],&Indices[0]);
#endif
  }
  Matrix->FillComplete();

  // ================================= //
  // print sparsity of original matrix //
  // ================================= //
 
  cout << "Sparsity, non-dropped matrix" << endl;
  Ifpack_PrintSparsity_Simple(*Matrix);

  // ====================================== //
  // create a new matrix, dropping by value //
  // ====================================== //
  //
  // drop all elements below 4.0. Only the upper-right element
  // is maintained, plus all the diagonals that are not
  // considering in dropping.
  Ifpack_DropFilter DropA(Matrix,4.0);
  assert (DropA.MaxNumEntries() == 2);

  cout << "Sparsity, dropping by value" << endl;
  Ifpack_PrintSparsity_Simple(DropA);

  // ========================================= //
  // create a new matrix, dropping by sparsity //
  // ========================================= //
  //
  // Mantain 2 off-diagonal elements.
  Ifpack_SparsityFilter SparsityA(Matrix,2);

  cout << "Sparsity, dropping by sparsity" << endl;
  Ifpack_PrintSparsity_Simple(SparsityA);
  assert (SparsityA.MaxNumEntries() == 3);

  // ======================================== //
  // create new matrices, dropping singletons //
  // ======================================== //
  //
  // If we apply this filter NumPoints - 1 times, 
  // we end up with a one-row matrix
  Ifpack_SingletonFilter Filter1(Matrix);
  Ifpack_SingletonFilter Filter2(Teuchos::rcp(&Filter1, false));
  Ifpack_SingletonFilter Filter3(Teuchos::rcp(&Filter2, false));
  Ifpack_SingletonFilter Filter4(Teuchos::rcp(&Filter3, false));

  cout << "Sparsity, dropping singletons 4 times" << endl;
  Ifpack_PrintSparsity_Simple(Filter4);
  assert (Filter4.NumMyRows() == 1);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return(EXIT_SUCCESS);
}

