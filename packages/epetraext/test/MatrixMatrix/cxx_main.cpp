//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
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
//@HEADER

//MatrixMatrix Test routine
#include <Epetra_ConfigDefs.h>
#include "EpetraExt_Version.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif

#include "Epetra_SerialComm.h"
#include "Epetra_Time.h"
#include "Epetra_Import.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "EpetraExt_MatrixMatrix.h"
#include "../epetra_test_err.h"

int check_rectangular_matrixmatrix_product_transAtransB(const Epetra_Comm& Comm, bool verbose);

int check_matrixmatrix_product(const Epetra_Comm& Comm, bool verbose);

int check_transpose_matrixmatrix_product(const Epetra_Comm& Comm, bool verbose);

int check_sparsedot();

int check_matrixmatrix_sum(const Epetra_Comm& Comm, bool verbose);

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

  if (verbose && Comm.MyPID()==0)
    cout << EpetraExt::EpetraExt_Version() << endl << endl;

  ierr = check_sparsedot();
  if (ierr != 0) {
    if (verbose && Comm.MyPID()==0) {
      cout << " sparsedot returned error code "<<ierr <<endl;
    }
  }

  if (Comm.NumProc() < 2) {
    ierr = check_rectangular_matrixmatrix_product_transAtransB(Comm, verbose);
    if (ierr != 0) {
      if (verbose && Comm.MyPID()==0) {
        cout << " rectangular matrix-matrix product returned code "<<ierr<<endl;
      }
    }
    else {
      if (verbose && Comm.MyPID()==0) {
        cout << " rectangular matrix-matrix product test passed"<<endl;
      }
    }
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

  ierr = check_transpose_matrixmatrix_product(Comm, verbose);
  if (ierr != 0) {
    if (verbose && Comm.MyPID()==0) {
      cout << " transpose matrix-matrix product returned code "<<ierr<<endl;
    }
  }
  else {
    if (verbose && Comm.MyPID()==0) {
      cout << " transpose matrix-matrix product test passed"<<endl;
    }
  }

  ierr = check_matrixmatrix_sum(Comm, verbose);
  if (ierr != 0) {
    if (verbose && Comm.MyPID()==0) {
      cout << " matrix-matrix sum returned code "<<ierr<<endl;
    }
  }
  else {
    if (verbose && Comm.MyPID()==0) {
      cout << " matrix-matrix sum test passed"<<endl;
    }
  }

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return returnierr;
}

int check_rectangular_matrixmatrix_product_transAtransB(const Epetra_Comm& Comm, bool verbose)
{
  int numProcs = Comm.NumProc();
  int myPID = Comm.MyPID();

  int numLocalRows = 5;
  int numGlobalRows = numLocalRows*numProcs;

  int numLocalCols = 4;
  int numGlobalCols = numLocalCols*numProcs;
  int indexBase = 0;

  Epetra_Map rowmap(-1, numLocalRows, indexBase, Comm);

  Epetra_Map colmap(-1, numLocalCols, indexBase, Comm);

  Epetra_CrsMatrix A(Copy, rowmap, colmap, 1);
  Epetra_CrsMatrix B(Copy, colmap, rowmap, 1);

  //For transAtransB, C will have A's row map and B's column map, which is
  //'colmap' in both cases (referring to the construction of A and B above).
  Epetra_CrsMatrix C(Copy, colmap, colmap, 1);

  int firstGlobalRow = numLocalRows*myPID;
  int i;

  //For this simple test, A and B are banded matrices filled with 1's, but they
  //are not symmetric. The bands are shifted to one side (up) of the diagonal.
  //The bands are of width 3 and include the diagonal and two positions above the
  //diagonal.

  numLocalRows = rowmap.NumMyElements();
  numLocalCols = colmap.NumMyElements();
  numGlobalRows = rowmap.NumGlobalElements();
  numGlobalCols = colmap.NumGlobalElements();

  for(i=0; i<numLocalRows; ++i) {
    int row = firstGlobalRow+i;

    int col = row;
    double value = 1.0;
    if (col < numGlobalCols) {
      A.InsertGlobalValues(row, 1, &value, &col);
    }

    if (col-1 < numGlobalCols && col-1 >= 0) {
      int colm1 = col-1;
      A.InsertGlobalValues(row, 1, &value, &colm1);
    }

    if (col < numGlobalCols-1) {
      int colp1 = col + 1;
      A.InsertGlobalValues(row, 1, &value, &colp1);
    }

    if (col < numGlobalCols-2) {
      int colp2 = col + 2;
      A.InsertGlobalValues(row, 1, &value, &colp2);
    }
  }

  numLocalCols = rowmap.NumMyElements();
  numLocalRows = colmap.NumMyElements();
  numGlobalCols = rowmap.NumGlobalElements();
  numGlobalRows = colmap.NumGlobalElements();

  for(i=0; i<numLocalRows; ++i) {
    int row = firstGlobalRow+i;

    int col = row;
    double value = 1.0;
    if (col < numGlobalCols) {
      B.InsertGlobalValues(row, 1, &value, &col);
    }

    if (col-1 < numGlobalCols && col-1 >= 0) {
      int colm1 = col-1;
      B.InsertGlobalValues(row, 1, &value, &colm1);
    }

    if (col < numGlobalCols-1) {
      int colp1 = col + 1;
      B.InsertGlobalValues(col, 1, &value, &colp1);
    }

    if (col < numGlobalCols-2) {
      int colp2 = col + 2;
      B.InsertGlobalValues(row, 1, &value, &colp2);
    }
  }

  EPETRA_CHK_ERR( A.FillComplete(rowmap, colmap) );
  EPETRA_CHK_ERR( B.FillComplete(colmap, rowmap) );

  if (verbose) {
     cout << "********** A **********"<<endl<<A<<endl;
  }

  if (verbose) {
     cout << "********** B **********"<<endl<<B<<endl;
  }

  EPETRA_CHK_ERR( EpetraExt::MatrixMatrix::Multiply(A, true, B, true, C) );

  if (verbose) {
    cout << "********** C = A*B **********"<<endl<<C<<endl;
  }

  return(0);
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

  //For this simple test, A and B are banded matrices filled with 1's, but they
  //are not symmetric. The bands are shifted to one side (up) of the diagonal.
  //The bands are of width 3 and include the diagonal and two positions above the
  //diagonal.

  for(i=0; i<numLocalRows; ++i) {
    int row = firstGlobalRow+i;

    int col = row;
    double value = 1.0;
    A.InsertGlobalValues(row, 1, &value, &col);
    B.InsertGlobalValues(row, 1, &value, &col);

    if (row < numGlobalRows-1) {
      int colp1 = col + 1;
      A.InsertGlobalValues(row, 1, &value, &colp1);
      B.InsertGlobalValues(row, 1, &value, &colp1);
    }

    if (row < numGlobalRows-2) {
      int colp2 = col + 2;
      A.InsertGlobalValues(row, 1, &value, &colp2);
      B.InsertGlobalValues(row, 1, &value, &colp2);
    }
  }

  EPETRA_CHK_ERR( A.FillComplete() );
  EPETRA_CHK_ERR( B.FillComplete() );

  EPETRA_CHK_ERR( EpetraExt::MatrixMatrix::Multiply(A, false, B, false, C) );

  if (verbose) {
//     cout << "********** A **********"<<endl<<A<<endl;
//     cout << "********** B **********"<<endl<<B<<endl;
    cout << "********** C = A*B **********"<<endl<<C<<endl;
  }

  return(0);
}

int check_transpose_matrixmatrix_product(const Epetra_Comm& Comm, bool verbose)
{
  int numProcs = Comm.NumProc();
  int myPID = Comm.MyPID();

  int numLocalRows = 4;
  int numGlobalRows = numLocalRows*numProcs;
  int indexBase = 0;

  Epetra_Map map(-1, numLocalRows, indexBase, Comm);

  Epetra_CrsMatrix A(Copy, map, 1);
  Epetra_CrsMatrix B(Copy, map, 1);
  Epetra_CrsMatrix C1(Copy, map, 1);
  Epetra_CrsMatrix C2(Copy, map, 1);
  Epetra_CrsMatrix C3(Copy, map, 1);

  int firstGlobalRow = numLocalRows*myPID;
  int i;

  //For this simple test, A and B are banded matrices filled with 1's, but they
  //are not symmetric. The bands are shifted to one side (up) of the diagonal.
  //The bands are of width 3 and include the diagonal and two positions above the
  //diagonal.

  for(i=0; i<numLocalRows; ++i) {
    int row = firstGlobalRow+i;

    int col = row;
    double value = 1.0;
    A.InsertGlobalValues(row, 1, &value, &col);
    B.InsertGlobalValues(row, 1, &value, &col);

    if (row < numGlobalRows-1) {
      int colp1 = col + 1;
      A.InsertGlobalValues(row, 1, &value, &colp1);
      B.InsertGlobalValues(row, 1, &value, &colp1);
    }

    if (row < numGlobalRows-2) {
      int colp2 = col + 2;
      A.InsertGlobalValues(row, 1, &value, &colp2);
      B.InsertGlobalValues(row, 1, &value, &colp2);
    }
  }

  EPETRA_CHK_ERR( A.FillComplete() );
  EPETRA_CHK_ERR( B.FillComplete() );

  EPETRA_CHK_ERR( EpetraExt::MatrixMatrix::Multiply(A, false, B, true, C1) );
  if (verbose) {
    cout << "C = A*B^T"<<endl;
    cout << C1 << endl;
    cout << endl<<"next, C = A^T*B"<<endl;
  }
  EPETRA_CHK_ERR( EpetraExt::MatrixMatrix::Multiply(A, true, B, false, C2) );
  if (verbose) {
    cout << C2 << endl;
    cout << endl << "next, C = A^T*B^T"<<endl;
  }
  EPETRA_CHK_ERR( EpetraExt::MatrixMatrix::Multiply(A, true, B, true, C3) );
  if (verbose) {
    cout << C3 << endl;
  }
  return(0);
}

int check_sparsedot()
{
  int u_len = 4;
  int* u_ind = new int[u_len];

  int v_len = 4;
  int* v_ind = new int[v_len];

  double* u = new double[u_len];
  double* v = new double[v_len];

  for(int i=0; i<u_len; ++i) {
    u[i] = 1.0*i;
    u_ind[i] = i;
  }

  for(int j=0; j<v_len; ++j) {
    v[j] = 1.0*j;
    v_ind[j] = j;
  }

  if (EpetraExt::sparsedot(u, u_ind, u_len,
			   v, v_ind, v_len) != 14.0) {
    return(-1);
  }

  u_ind[0] = 3;
  u_ind[1] = 6;
  u_ind[2] = 7;
  u_ind[3] = 9;

  v_ind[0] = 1;
  v_ind[1] = 3;
  v_ind[2] = 6;
  v_ind[3] = 10;

  if (EpetraExt::sparsedot(u, u_ind, u_len,
			   v, v_ind, v_len) != 2.0) {
    return(-1);
  }

  delete [] u_ind;
  delete [] v_ind;
  delete [] u;
  delete [] v;

  return(0);
}

int check_matrixmatrix_sum(const Epetra_Comm& Comm, bool verbose)
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

  //For this simple test, A and B are banded matrices filled with 1's, but they
  //are not symmetric. The bands are shifted to one side (up) of the diagonal.
  //The bands are of width 3 and include the diagonal and two positions above the
  //diagonal.

  for(i=0; i<numLocalRows; ++i) {
    int row = firstGlobalRow+i;

    int col = row;
    double value = 1.0;
    A.InsertGlobalValues(row, 1, &value, &col);
    B.InsertGlobalValues(row, 1, &value, &col);

    if (row < numGlobalRows-1) {
      int colp1 = col + 1;
      A.InsertGlobalValues(row, 1, &value, &colp1);
      B.InsertGlobalValues(row, 1, &value, &colp1);
    }

    if (row < numGlobalRows-2) {
      int colp2 = col + 2;
      A.InsertGlobalValues(row, 1, &value, &colp2);
      B.InsertGlobalValues(row, 1, &value, &colp2);
    }
  }

  EPETRA_CHK_ERR( A.FillComplete() );
  EPETRA_CHK_ERR( B.FillComplete() );

  EPETRA_CHK_ERR( EpetraExt::MatrixMatrix::Add(A, true, 1.0, C, 0.0) );
  EPETRA_CHK_ERR( EpetraExt::MatrixMatrix::Add(B, false, 1.0, C, 1.0) );

  if (verbose) {
    cout << "********** C = (B + (A^T + C)) **********"<<C<<endl;
  }

  return(0);
}

