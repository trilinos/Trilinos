#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#ifdef PETRA_MPI
#include "mpi.h"
#endif
#include "Trilinos_Util.h"
#ifndef __cplusplus
#define __cplusplus
#endif
#include "Petra_Comm.h"
#include "Petra_Time.h"
#include "Petra_Map.h"
#include "Petra_BlockMap.h"
#include "Petra_RDP_MultiVector.h"
#include "Petra_RDP_Vector.h"
#include "Petra_RDP_CRS_Matrix.h"
#include "Petra_CRS_Graph.h"

int main(int argc, char *argv[])
{
  
  int ierr = 0, i, j;
  bool debug = false;

#ifdef PETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;

#endif
  bool verbose = true;
  /*
  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;
  */


#ifdef PETRA_MPI
  Petra_Comm & Comm = *new Petra_Comm( MPI_COMM_WORLD );
#else
  Petra_Comm & Comm = *new Petra_Comm();
#endif

 /* 
  char tmp;
  if (rank==0) cout << "Press any key to continue..."<< endl;
  if (rank==0) cin >> tmp;
  Comm.Barrier();
 */ 
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();
  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
              << " is alive."<<endl;

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


  Petra_Map& Map = *new Petra_Map(NumGlobalEquations, NumMyEquations, 
			MyGlobalEquations, 0, Comm);
 
  Petra_Time & timer = *new Petra_Time(Comm);

  double start = timer.ElapsedTime();
  Petra_RDP_CRS_Matrix& A = *new Petra_RDP_CRS_Matrix(Copy, Map, NumNz);
  
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
  
  assert(A.TransformToLocal()==0);
  
  if (verbose) cout << "\nTime to construct A                = " << timer.ElapsedTime() - start << endl;
  double * xexactt = xexact;
  Petra_RDP_Vector& xx = *new Petra_RDP_Vector(Copy, Map, xexactt);

  double * bt = b;
  Petra_RDP_Vector& bb = *new Petra_RDP_Vector(Copy, Map, bt);
  Petra_RDP_Vector& bcomp = *new Petra_RDP_Vector(Map);

  // Sanity check
  assert(A.Multiply(false, xx, bcomp)==0);
  Petra_RDP_Vector& resid = *new Petra_RDP_Vector(Map); 
 
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
  Petra_CRS_Graph AG = A.Graph(); // Get matrix graph

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

  const Petra_Map & TransMap = A.ImportMap();

  Petra_RDP_CRS_Matrix & TempTransA1 = *new Petra_RDP_CRS_Matrix(View, TransMap, TransNumNz);
  int * TransMyGlobalEquations = new int[NumMyCols];
  TransMap.MyGlobalElements(TransMyGlobalEquations);
  
  /* Add  rows one-at-a-time */

  for (i=0; i<NumMyCols; i++)
    {
     assert(TempTransA1.InsertGlobalValues(TransMyGlobalEquations[i], 
				      TransNumNz[i], TransValues[i], TransIndices[i])==0);
    }
 
  // Note: The following call to TransformToLocal is currently necessary because
  //       some global constants that are needed by the Export () are computed in this routine
  assert(TempTransA1.TransformToLocal()==0);

  // Now that transpose matrix with shared rows is entered, create a new matrix that will
  // get the transpose with uniquely owned rows (using the same row distribution as A).

  Petra_RDP_CRS_Matrix& TransA1 = *new Petra_RDP_CRS_Matrix(Copy, Map,0);

  // Create an Export object that will move TempTransA around

  Petra_Export & Export = *new Petra_Export(TransMap, Map);

  assert(TransA1.Export(TempTransA1, Export, Add)==0);
  
  assert(TransA1.TransformToLocal()==0);


  if (verbose) cout << "\nTime to construct TransA1          = " << timer.ElapsedTime() - start << endl;

  // Now compute b = A' * x using the transpose option on Multiply and using 
  // created transpose matrix

  Petra_RDP_Vector& x = *new Petra_RDP_Vector(Map);
  x.Random(); // Fill with random numbers

  Petra_RDP_Vector& b1 = *new Petra_RDP_Vector(Map);
  assert(A.Multiply(true, x, b1)==0);
  Petra_RDP_Vector& b2 = *new Petra_RDP_Vector(Map);
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

  // const Petra_Map & TransMap = A.ImportMap();

  start = timer.ElapsedTime();

  Petra_RDP_CRS_Matrix & TempTransA2 = *new Petra_RDP_CRS_Matrix(Copy, TransMap, 0);
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


  
  // Note: The following call to TransformToLocal is currently necessary because
  //       some global constants that are needed by the Export () are computed in this routine
  assert(TempTransA2.TransformToLocal()==0);

  if (verbose) cout << "\nTime to construct TransA2          = " << timer.ElapsedTime() - start << endl;

  // Now that transpose matrix with shared rows is entered, create a new matrix that will
  // get the transpose with uniquely owned rows (using the same row distribution as A).

  Petra_RDP_CRS_Matrix& TransA2 = *new Petra_RDP_CRS_Matrix(Copy, Map,0);

  // Create an Export object that will move TempTransA around

  // Petra_Export & Export = *new Petra_Export(TransMap, Map); // Export already built

  assert(TransA2.Export(TempTransA2, Export, Add)==0);
  
  assert(TransA2.TransformToLocal()==0);

  // Now compute b = A' * x using the transpose option on Multiply and using 
  // created transpose matrix

  // Petra_RDP_Vector& x = *new Petra_RDP_Vector(Map);
  // x.Random(); // Fill with random numbers

  // Petra_RDP_Vector& b1 = *new Petra_RDP_Vector(Map);
  assert(A.Multiply(true, x, b1)==0);
  // Petra_RDP_Vector& b2 = *new Petra_RDP_Vector(Map);
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

  delete &bcomp;
  delete &bb;
  delete &xx;
  delete &b1;
  delete &b2;
  delete &resid;
  delete &x;
  delete &A;
  delete &TransA1;
  delete &TransA2;
  delete &TempTransA1;
  delete &TempTransA2;

  delete &Export;

  delete &Map;
  delete &Comm;
				       
#ifdef PETRA_MPI
  MPI_Finalize() ;
#endif

return 0 ;
}
