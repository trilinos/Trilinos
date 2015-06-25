/*@HEADER
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
*/

#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Ifpack_ILU.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_AdditiveSchwarz.h"

template<class T>
void Test(const std::string what, Epetra_RowMatrix& A)
{
  T Prec(&A);

  bool UseTranspose = true;

  IFPACK_CHK_ERRV(Prec.Initialize());
  IFPACK_CHK_ERRV(Prec.Compute());
  IFPACK_CHK_ERRV(Prec.SetUseTranspose(UseTranspose));

  Epetra_MultiVector LHS_exact(A.OperatorDomainMap(), 2);
  Epetra_MultiVector LHS(A.OperatorDomainMap(), 2);
  Epetra_MultiVector RHS(A.OperatorRangeMap(), 2);

  LHS_exact.Random(); LHS.PutScalar(0.0);

  A.Multiply(UseTranspose, LHS_exact, RHS);

  Prec.ApplyInverse(RHS, LHS);

  LHS.Update(1.0, LHS_exact, -1.0);
  double norm[2];

  LHS.Norm2(norm);
  norm[0] += norm[1];

  if (norm[0] > 1e-5)
  {
    cout << what << ": Test failed: norm = " << norm[0] << endl;
    exit(EXIT_FAILURE);
  }

  cout << what << ": Test passed: norm = " << norm[0] << endl;
}

// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[]) 
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  if (Comm.NumProc() != 1)
  {
    cerr << "To be run with one processor only" << endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_SUCCESS);
  }

  Epetra_Map Map(8, 0, Comm);

  Epetra_CrsMatrix A(Copy, Map, 0);

  // for this matrix the incomplete factorization
  // is the exact one, so ILU and ILUT must be exact solvers.
  for (int row = 0; row < 8; ++row)
  {
    double value = 2.0 + row;
    A.InsertGlobalValues(row, 1, &value, &row);
    if (row)
    {
      int col = row - 1;
      value = 1.0 + row;
      A.InsertGlobalValues(row, 1, &value, &col);
    }
#if 0
    if (row != Map.NumGlobalElements() - 1)
    {
      int col = row + 1;
      value = 0.0;
      A.InsertGlobalValues(row, 1, &value, &col);
    }
#endif
  }

  A.FillComplete();

  Test<Ifpack_ILU>("Ifpack_ILU", A);
  Test<Ifpack_ILUT>("Ifpack_ILUT", A);
  Test<Ifpack_AdditiveSchwarz<Ifpack_ILU> >("AS, Ifpack_ILU", A);
  Test<Ifpack_AdditiveSchwarz<Ifpack_ILUT> >("AS, Ifpack_ILUT", A);

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return(EXIT_SUCCESS);
}
