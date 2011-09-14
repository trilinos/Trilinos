//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
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


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "Epetra_Comm.h"
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Export.h"
#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Trilinos_Util.h"
#ifndef __cplusplus
#define __cplusplus
#endif

int main(int argc, char *argv[])
{
  
  int ierr = 0, i, j;
  bool debug = false;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;
  Epetra_SerialComm Comm;

#endif
  bool verbose = true;
  /*
  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;
  */


 
  char tmp;
  if (rank==0) cout << "Press any key to continue..."<< endl;
  if (rank==0) cin >> tmp;
  Comm.Barrier();

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();
  if (verbose) cout << Comm << endl;

  bool verbose1 = verbose;

  // Redefine verbose to only print on PE 0
  if (verbose && rank!=0) verbose = false;


  if(argc != 2) cout << "Error: Enter name of data file on command line." << endl; 

  /* Read matrix file and distribute among processors.  
     Returns with this processor's set of rows */ 

  int NumGlobalEquations, n_nonzeros, *bindx;
  double * val, * xguess, * b, * xexact = 0;

  Trilinos_Util_read_hb(argv[1], Comm.MyPID(), &NumGlobalEquations, &n_nonzeros,
			&val,  &bindx, &xguess, &b, &xexact);

  int NumMyEquations, * MyGlobalEquations;

  Trilinos_Util_distrib_msr_matrix(Comm, &NumGlobalEquations, &n_nonzeros, &NumMyEquations,
		  &MyGlobalEquations, &val, &bindx, &xguess, &b, &xexact);


  /* Make NumNz - number of entries in each row */

  int * NumNz = new int[NumMyEquations];
  for (i=0; i<NumMyEquations; i++) NumNz[i] = bindx[i+1] - bindx[i] + 1;


  Epetra_Map Map(NumGlobalEquations, NumMyEquations, 
			MyGlobalEquations, 0, Comm);
 
  Epetra_Time timer(Comm);

  double start = timer.ElapsedTime();
  Epetra_CrsMatrix A(Copy, Map, NumNz);
  
  /* Add  rows one-at-a-time */

  int NumIndices;
  int * Indices;
  double * Values;
  for (i=0; i<NumMyEquations; i++){
      Values = val + bindx[i];
      Indices = bindx + bindx[i];
      NumIndices = bindx[i+1] - bindx[i];
     assert(A.InsertGlobalValues(MyGlobalEquations[i], NumIndices, Values, Indices)==0);
     assert(A.InsertGlobalValues(MyGlobalEquations[i], 1, val+i, MyGlobalEquations+i)==0);
    }
  
  assert(A.FillComplete()==0);
  
  if (verbose) cout << "\nTime to construct A                = " << timer.ElapsedTime() - start << endl;
  double * xexactt = xexact;
  Epetra_Vector xx(Copy, Map, xexactt);

  double * bt = b;
  Epetra_Vector bb(Copy, Map, bt);
  Epetra_Vector bcomp(Map);

  // Sanity check
  assert(A.Multiply(false, xx, bcomp)==0);
  Epetra_Vector resid(Map); 
 
  assert(resid.Update(1.0, bb, -1.0, bcomp, 0.0)==0);

  double residual;
  assert(resid.Norm2(&residual)==0);
  if (Comm.MyPID()==0) cout << "Sanity check: Norm of b - Ax for known x and b = " << residual << endl;

  // Way 1: This approach is probably most time-efficient, but is a little more complex because
  //        we explicitly pre-compute the local transpose.  It should not use anymore memory
  //        than Way 2 since we make a view of the pre-export local transpose, which we
  //        cannot do in Way 2.

  // Extract newly constructed matrix one row at a time
  // 1) First get the local indices to count how many nonzeros will be in the 
  //    transpose graph on each processor

  start = timer.ElapsedTime();
  Epetra_CrsGraph AG = A.Graph(); // Get matrix graph

  int NumMyCols = AG.NumMyCols();
  int * TransNumNz = new int[NumMyCols];
  for (i=0;i<NumMyCols; i++) TransNumNz[i] = 0;
  for (i=0; i<NumMyEquations; i++) {
    assert(AG.ExtractMyRowView(i, NumIndices, Indices)==0); // Get view of ith row
    for (j=0; j<NumIndices; j++) ++TransNumNz[Indices[j]];
  }

  int ** TransIndices = new int*[NumMyCols];
  double ** TransValues = new double*[NumMyCols];

  for(i=0; i<NumMyCols; i++) {
    NumIndices = TransNumNz[i];
    if (NumIndices>0) {
      TransIndices[i] = new int[NumIndices];
      TransValues[i] = new double[NumIndices];
    }
  }

  // Now copy values and global indices into newly create transpose storage

  for (i=0;i<NumMyCols; i++) TransNumNz[i] = 0; // Reset transpose NumNz counter
  for (i=0; i<NumMyEquations; i++) {
    assert(A.ExtractMyRowView(i, NumIndices, Values, Indices)==0);
    int ii = A.GRID(i);
    for (j=0; j<NumIndices; j++) {
      int TransRow = Indices[j];
      int loc = TransNumNz[TransRow];
      TransIndices[TransRow][loc] = ii;
      TransValues[TransRow][loc] = Values[j];
      ++TransNumNz[TransRow]; // increment counter into current transpose row
    }
  }

  //  Build Transpose matrix with some rows being shared across processors.
  //  We will use a view here since the matrix will not be used for anything else

  const Epetra_Map & TransMap = A.ImportMap();

  Epetra_CrsMatrix TempTransA1(View, TransMap, TransNumNz);
  int * TransMyGlobalEquations = new int[NumMyCols];
  TransMap.MyGlobalElements(TransMyGlobalEquations);
  
  /* Add  rows one-at-a-time */

  for (i=0; i<NumMyCols; i++)
    {
     assert(TempTransA1.InsertGlobalValues(TransMyGlobalEquations[i], 
				      TransNumNz[i], TransValues[i], TransIndices[i])==0);
    }
 
  // Note: The following call to FillComplete is currently necessary because
  //       some global constants that are needed by the Export () are computed in this routine
  assert(TempTransA1.FillComplete()==0);

  // Now that transpose matrix with shared rows is entered, create a new matrix that will
  // get the transpose with uniquely owned rows (using the same row distribution as A).

  Epetra_CrsMatrix TransA1(Copy, Map,0);

  // Create an Export object that will move TempTransA around

  Epetra_Export Export(TransMap, Map);

  assert(TransA1.Export(TempTransA1, Export, Add)==0);
  
  assert(TransA1.FillComplete()==0);


  if (verbose) cout << "\nTime to construct TransA1          = " << timer.ElapsedTime() - start << endl;

  // Now compute b = A' * x using the transpose option on Multiply and using 
  // created transpose matrix

  Epetra_Vector x(Map);
  x.Random(); // Fill with random numbers

  Epetra_Vector b1(Map);
  assert(A.Multiply(true, x, b1)==0);
  Epetra_Vector b2(Map);
  assert(TransA1.Multiply(false, x, b2)==0);
 
  assert(resid.Update(1.0, b1, -1.0, b2, 0.0)==0);

  assert(b1.Norm2(&residual)==0);
  if (verbose) cout << "Norm of RHS using Trans = true with A           = " << residual << endl;
  assert(b2.Norm2(&residual)==0);
  if (verbose) cout << "Norm of RHS using Trans = false with TransA1    = " << residual << endl;
  assert(resid.Norm2(&residual)==0);
  if (verbose) cout << "Difference between using A and TransA1          = " << residual << endl;


  // Way 2: This approach is probably the simplest to code, but is probably slower.
  //        We compute the transpose by dumping entries one-at-a-time.  

  // Extract newly constructed matrix one entry at a time and
  //  build Transpose matrix with some rows being shared across processors.

  // const Epetra_Map & TransMap = A.ImportMap();

  start = timer.ElapsedTime();

  Epetra_CrsMatrix TempTransA2(Copy, TransMap, 0);
  TransMap.MyGlobalElements(TransMyGlobalEquations);

  for (int LocalRow=0; LocalRow<NumMyEquations; LocalRow++) {
    assert(A.ExtractMyRowView(LocalRow, NumIndices, Values, Indices)==0);
    int TransGlobalCol = A.GRID(LocalRow);
    for (j=0; j<NumIndices; j++) {
      int TransGlobalRow = A.GCID(Indices[j]);
      double TransValue = Values[j];
      assert(TempTransA2.InsertGlobalValues(TransGlobalRow, 1, &TransValue, &TransGlobalCol)>=0);
    }
  }


  
  // Note: The following call to FillComplete is currently necessary because
  //       some global constants that are needed by the Export () are computed in this routine
  assert(TempTransA2.FillComplete()==0);

  if (verbose) cout << "\nTime to construct TransA2          = " << timer.ElapsedTime() - start << endl;

  // Now that transpose matrix with shared rows is entered, create a new matrix that will
  // get the transpose with uniquely owned rows (using the same row distribution as A).

  Epetra_CrsMatrix TransA2(Copy, Map,0);

  // Create an Export object that will move TempTransA around

  // Epetra_Export Export(TransMap, Map); // Export already built

  assert(TransA2.Export(TempTransA2, Export, Add)==0);
  
  assert(TransA2.FillComplete()==0);

  // Now compute b = A' * x using the transpose option on Multiply and using 
  // created transpose matrix

  // Epetra_Vector x(Map);
  // x.Random(); // Fill with random numbers

  // Epetra_Vector b1(Map);
  assert(A.Multiply(true, x, b1)==0);
  // Epetra_Vector b2(Map);
  assert(TransA2.Multiply(false, x, b2)==0);
 
  assert(resid.Update(1.0, b1, -1.0, b2, 0.0)==0);

  assert(b1.Norm2(&residual)==0);
  if (verbose) cout << "Norm of RHS using Trans = true with A           = " << residual << endl;
  assert(b2.Norm2(&residual)==0);
  if (verbose) cout << "Norm of RHS using Trans = false with TransA2    = " << residual << endl;
  assert(resid.Norm2(&residual)==0);
  if (verbose) cout << "Difference between using A and TransA2          = " << residual << endl;


  // The free's are needed because the HB utility routines still C-style calloc calls
  free ((void *) xguess);
  free ((void *) b);
  free ((void *) xexact);
  free ((void *) val);
  free ((void *) bindx);
  free ((void *) MyGlobalEquations);

  delete [] TransMyGlobalEquations;

  for(i=0; i<NumMyCols; i++) {
    NumIndices = TransNumNz[i];
    if (NumIndices>0) {
      delete [] TransIndices[i];
      delete [] TransValues[i];
    }
  }
  delete [] TransIndices;
  delete [] TransValues;

  delete [] NumNz;
  delete [] TransNumNz;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

return 0 ;
}
