//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

//MatrixMatrix Test routine
#include <Epetra_ConfigDefs.h>

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_Time.h"
#include "Epetra_Import.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "EpetraExt_MatrixMatrix.h"
#include "../epetra_test_err.h"

int check_matrixmatrix_product(const Epetra_Comm& Comm, bool verbose);

int main(int argc, char *argv[]) {

  int ierr=0, returnierr=0;

  bool verbose = false;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

#else
  Epetra_SerialComm Comm;
#endif

  // Check if we should print results to standard out
  if (argc>1) {
    if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;
  }

  //Make sure the value of verbose is consistent across processors.
  int verbose_int = verbose ? 1 : 0;
  Comm.Broadcast(&verbose_int, 1, 0);
  verbose = verbose_int==1 ? true : false;

  if (!verbose) {
    Comm.SetTracebackMode(0); // This should shut down error traceback reporting
  }

  ierr = check_matrixmatrix_product(Comm, verbose);
  if (ierr != 0) {
    if (verbose && Comm.MyPID()==0) {
      cout << " matrix-matrix product returned code "<<ierr<<endl;
    }
  }
  else {
    if (verbose && Comm.MyPID()==0) {
      cout << " matrix-matrix product test passed"<<endl;
    }
  }

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return returnierr;
}

int check_matrixmatrix_product(const Epetra_Comm& Comm, bool verbose)
{
  int numProcs = Comm.NumProc();
  int myPID = Comm.MyPID();

  int numLocalRows = 4;
  int numGlobalRows = numLocalRows*numProcs;
  int indexBase = 0;

  Epetra_Map map(-1, numLocalRows, indexBase, Comm);

  Epetra_CrsMatrix A(Copy, map, 1);
  Epetra_CrsMatrix B(Copy, map, 1);
  Epetra_CrsMatrix C(Copy, map, 1);

  int firstGlobalRow = numLocalRows*myPID;
  int i;

  //For this simple test, A and B are tri-diagonal matrices filled with 1's.

  for(i=0; i<numLocalRows; ++i) {
    int row = firstGlobalRow+i;

    int col = row;
    double value = 1.0;
    //double bvalue = value*(row+1);
    A.InsertGlobalValues(row, 1, &value, &col);
    B.InsertGlobalValues(row, 1, &value, &col);

    if (row > 0) {
      int colm1 = col - 1;
      //bvalue -= 0.1;
      A.InsertGlobalValues(row, 1, &value, &colm1);
      B.InsertGlobalValues(row, 1, &value, &colm1);
      //bvalue += 0.1;
    }

    if (row < numGlobalRows-1) {
      int colp1 = col + 1;
      //bvalue += 0.1;
      A.InsertGlobalValues(row, 1, &value, &colp1);
      B.InsertGlobalValues(row, 1, &value, &colp1);
    }
  }

  EPETRA_CHK_ERR( A.FillComplete() );
  EPETRA_CHK_ERR( B.FillComplete() );

  EPETRA_CHK_ERR( EpetraExt::MatrixMatrix::Product(A, B, C) );

  //For these simple operands, the result C should have a simple form which
  //we can check as follows.

  int len = 20;
  int numIndices;
  int* indices = new int[len];
  double* values = new double[len];

  for(i=0; i<numLocalRows; ++i) {
    int row = firstGlobalRow+i;

    EPETRA_CHK_ERR( C.ExtractGlobalRowCopy(row, len, numIndices,
					   values, indices) );

    double sum = 0;
    for(int j=0; j<numIndices; ++j) {
      sum += values[j];
    }

    if (row == 0 || row == numGlobalRows-1) {
      if (numIndices != 3) {
	return(-2);
      }
      if (sum != 5.0) {
	return(-3);
      }
    }
    else if (row == 1 || row == numGlobalRows-2) {
      if (numIndices != 4) {
	cout << "row "<<row<<", numIndices ("<<numIndices<<") should be "<<4<<endl;
	return(-4);
      }
      if (sum != 8.0) {
 	return(-5);
      }
    }
    else {
      if (numIndices != 5) {
	cout << "row "<<row<<", numIndices ("<<numIndices<<") should be "<<5<<endl;
	return(-6);
      }
      if (sum != 9.0) {
 	return(-7);
      }
    }

  }

//  if (verbose) {
//     cout << "********** A **********"<<endl<<A<<endl;
//     cout << "********** B **********"<<endl<<B<<endl;
//     cout << "********** C **********"<<endl<<C<<endl;
//  }

  delete [] indices;
  delete [] values;

  return(0);
}
