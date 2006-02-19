// @HEADER
// ***********************************************************************
// 
//            Trilinos: An Object-Oriented Solver Framework
//                 Copyright (2001) Sandia Corporation
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
void Test(const string what, Epetra_RowMatrix& A)
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
