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
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Ifpack_Reordering.h"
#include "Ifpack_RCMReordering.h"
#include "Ifpack_ReorderFilter.h"
#include "Ifpack_Utils.h"

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
  
  Epetra_Map Map(-1,NumPoints,0,Comm);
  
  vector<int> Indices(NumPoints);
  vector<double> Values(NumPoints);

  Epetra_CrsMatrix A(Copy,Map,0);
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

    A.InsertGlobalValues(i, NumEntries, &Values[0], &Indices[0]);
  }

  A.FillComplete();

  // print the sparsity to file, postscript format
  ////Ifpack_PrintSparsity(A,"OrigA.ps");

  // create the reordering...
  Ifpack_RCMReordering Reorder;
  // and compute is on A
  IFPACK_CHK_ERR(Reorder.Compute(A));

  // cout information
  cout << Reorder;

  // create a reordered matrix
  Ifpack_ReorderFilter ReordA(&A,&Reorder);

  // print the sparsity to file, postscript format
  ////Ifpack_PrintSparsity(ReordA,"ReordA.ps");

#ifdef HAVE_MPI
  MPI_Finalize(); 
#endif
  return(EXIT_SUCCESS);

}
