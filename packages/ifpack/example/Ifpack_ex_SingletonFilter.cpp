// @HEADER
// ***********************************************************************
// 
//                IFPACK
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

#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Ifpack_SingletonFilter.h"
#if defined(HAVE_IFPACK_AMESOS)
#include "Ifpack_Amesos.h"
#endif

//==============================================================================
void PrintSparsity(Epetra_RowMatrix& A)
{
  int MaxEntries = A.MaxNumEntries();
  vector<int> Indices(MaxEntries);
  vector<double> Values(MaxEntries);
  vector<bool> FullRow(A.NumMyRows());

  cout << "NumMyRows() = " << A.NumMyRows() << endl;

  cout << "+";
  for (int j = 0 ; j < A.NumMyRows() ; ++j)
    cout << '-';
  cout << "+" << endl;

  for (int i = 0 ; i < A.NumMyRows() ; ++i) {

    int Length;
    A.ExtractMyRowCopy(i,MaxEntries,Length,
		       &Values[0], &Indices[0]);

    for (int j = 0 ; j < A.NumMyRows() ; ++j)
      FullRow[j] = false;
    
    for (int j = 0 ; j < Length ; ++j) {
      FullRow[Indices[j]] = true;
    }

    cout << "|";
    for (int j = 0 ; j < A.NumMyRows() ; ++j) {
      if (FullRow[j])
	cout << 'x';
      else
	cout << ' ';
    }
    cout << "|" << endl;
  }

  cout << "+";
  for (int j = 0 ; j < A.NumMyRows() ; ++j)
    cout << '-';
  cout << "+" << endl << endl;

}

//==============================================================================
int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

#ifdef HAVE_IFPACK_TEUCHOS

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

  int NumPoints = 16;
  
  Epetra_Map Map(-1,NumPoints,0,Comm);
  
  vector<int> Indices(NumPoints);
  vector<double> Values(NumPoints);

  Epetra_CrsMatrix A(Copy,Map,0);
  for (int i = 0 ; i < NumPoints ; ++i) {
    
    if ((i != 0) && (i != NumPoints - 1)) {
      Indices[0] = i - 1;
      Values[0] = -1.0;
      Indices[1] = i + 1;
      Values[1] = -1.0;
      A.InsertGlobalValues(i, 2, &Values[0], &Indices[0]);
    }
    // now diagonal
    Indices[0] = i;
    Values[0] = 6.0;
    A.InsertGlobalValues(i, 1, &Values[0], &Indices[0]);
  }

  A.FillComplete();

  cout << "Print the sparsity of original matrix" << endl;
  PrintSparsity(A);

  // create a wrapper that eliminates the singletons
  Ifpack_SingletonFilter FilterA(&A);

  cout << "Print the sparsity of singleton-filtered matrix" << endl;
  PrintSparsity(FilterA);

#ifdef PRINT_MATRIX
  Indices.resize(FilterA.MaxNumEntries());
  Values.resize(FilterA.MaxNumEntries());
  for (int i = 0 ; i < FilterA.NumMyRows() ; ++i) {

    int Length;
    FilterA.ExtractMyRowCopy(i,FilterA.MaxNumEntries(),Length,
		       &Values[0], &Indices[0]);
    for (int j = 0 ; j < Length ; ++j) {
      cout << "A(" << i << "," << Indices[j] << ") = " << Values[j] << endl;
    }
  }
#endif

#ifdef HAVE_IFPACK_AMESOS
  // attempt the solve Ax = b using the filtered matrix
  Epetra_Vector LHSexact(Map);
  for (int i = 0 ; i < NumPoints ; ++i)
    LHSexact[i] = 5.0 + 1.0 * i;
  Epetra_Vector LHS(Map);
  Epetra_Vector RHS(Map);
  A.Multiply(false,LHSexact,RHS);
  LHS.PutScalar(0.0);

  // use Amesos solver (through IFPACK)
  Ifpack_Amesos Solver(&FilterA);
  IFPACK_CHK_ERR(Solver.Initialize());
  IFPACK_CHK_ERR(Solver.Compute());

  Epetra_Vector ReducedRHS(FilterA.Map());
  Epetra_Vector ReducedLHS(FilterA.Map());
  FilterA.SolveSingletons(RHS,LHS);

  FilterA.CreateReducedRHS(LHS,RHS,ReducedRHS);
  IFPACK_CHK_ERR(Solver.ApplyInverse(ReducedRHS,ReducedLHS));

  FilterA.UpdateLHS(ReducedLHS,LHS);

  LHSexact.Update(1.0,LHS,-1.0);
  vector<double> Norm2(1);

  LHSexact.Norm2(&Norm2[0]);
  for (int i = 0 ; i < 1 ; ++i) {
    cout << Norm2[i] << endl;
    if (Norm2[i] > 1e-5) {
      cerr << "TEST FAILED!" << endl;
      exit(EXIT_FAILURE);
    }
  }
#endif

  cout << "TEST PASSED!" << endl;
#else
  puts("Please configure ifpack with --enable-teuchos to run this example");
#endif

#ifdef HAVE_MPI
  MPI_Finalize(); 
#endif
  return(EXIT_SUCCESS);

}
