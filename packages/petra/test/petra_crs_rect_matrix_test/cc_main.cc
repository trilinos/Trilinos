#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#ifdef PETRA_MPI
#include "mpi.h"
#endif
#ifndef __cplusplus
#define __cplusplus
#endif
#include "Petra_Comm.h"
#include "Petra_Map.h"
#include "Petra_Time.h"
#include "Petra_RDP_Vector.h"
#include "Petra_RDP_CRS_Matrix.h"

 
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

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;



#ifdef PETRA_MPI
  Petra_Comm & Comm = *new Petra_Comm( MPI_COMM_WORLD );
#else
  Petra_Comm & Comm = *new Petra_Comm();
#endif


  //  char tmp;
  //  if (rank==0) cout << "Press any key to continue..."<< endl;
  //  if (rank==0) cin >> tmp;
  //  Comm.Barrier();

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();
  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
              << " is alive."<<endl;

  bool verbose1 = verbose;

  // Redefine verbose to only print on PE 0
  if (verbose && rank!=0) verbose = false;

  int NumMyEquations = 10000;

  int NumGlobalEquations = NumMyEquations*NumProc;
  int NumGlobalVariables = 2 * NumGlobalEquations+1;

  // Construct a Map that puts approximately the same Number of equations on each processor

  Petra_Map& RowMap = *new Petra_Map(NumGlobalEquations, 0, Comm);
  Petra_Map& XMap = *new Petra_Map(NumGlobalVariables, 0, Comm);
  Petra_Map& YMap = RowMap;
  
  int NumMyY = NumMyEquations;
  int IndexBase = 0;
  bool DistributedGlobal = (NumGlobalEquations>NumMyEquations);

  // Get update list and number of local equations from newly created Map
  int * MyGlobalElements = new int[RowMap.NumMyElements()];
  RowMap.MyGlobalElements(MyGlobalElements);

  // Get update list and number of local equations from newly created XMap
  int * XGlobalElements = new int[XMap.NumMyElements()];
  XMap.MyGlobalElements(XGlobalElements);

  // Get update list and number of local variables from newly created YMap
  int * YGlobalElements = new int[YMap.NumMyElements()];
  YMap.MyGlobalElements(YGlobalElements);

  // We need vectors to compute:
  // X = A^T*Y
  // AATY = A*A^T*Y = A*X
  //  and 
  // BY = B*Y

  Petra_RDP_Vector & Y = *new Petra_RDP_Vector(YMap);
  Petra_RDP_Vector & X = *new Petra_RDP_Vector(XMap);
  Petra_RDP_Vector & AATY = *new Petra_RDP_Vector(YMap);
  Petra_RDP_Vector & BY = *new Petra_RDP_Vector(YMap);


  // Fill Y Vector
  Y.Random();

  // Create a Petra_Matrix with the values of A
  // A is a simple 1D weighted average operator that mimics a restriction operator
  // that might be found in a multigrid code.

  Petra_RDP_CRS_Matrix& A = *new Petra_RDP_CRS_Matrix(Copy, RowMap, 3);
  
  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1


  double *Values = new double[3];
  int *Indices = new int[3];
  int NumEntries;
  Values[0] = 0.25;
  Values[1] = 0.5;
  Values[2] = 0.25;
  for (i=0; i<NumMyEquations; i++)
    {
      Indices[0] = 2*MyGlobalElements[i];
      Indices[1] = 2*MyGlobalElements[i]+1;
      Indices[2] = 2*MyGlobalElements[i]+2;
      NumEntries = 3;
     assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices)==0);
    }
  
  // Finish up
  assert(A.TransformToLocal(&XMap, &YMap)==0);


  if (NumGlobalEquations<20) cout << "\n\n Matrix A = " << A << endl;


  // Create a Petra_Matrix containing B = A*A^T.
  // This matrix will be a square tridiagonal matrix.  We will use it to compare the results
  // of A*(A^T*X) using two methods: (1) with two calls to Multiply using A^T and then A and
  // (2) using B directly.

  Petra_RDP_CRS_Matrix& B = *new Petra_RDP_CRS_Matrix(Copy, RowMap, 3);

  Values[0] = 1.0/16.0;
  Values[1] = 3.0/8.0;
  Values[2] = 1.0/16.0;
  int Valstart;
  for (i=0; i<NumMyEquations; i++)
    {
      if (MyGlobalElements[i] == 0) {
      Indices[0] = MyGlobalElements[i];
      Indices[1] = MyGlobalElements[i]+1;
      NumEntries = 2;
      Valstart = 1;
      }
      else {
	Indices[0] = MyGlobalElements[i]-1;
	Indices[1] = MyGlobalElements[i];
	Indices[2] = MyGlobalElements[i]+1;
	NumEntries = 3;
	Valstart = 0;
      }
      if (MyGlobalElements[i] == NumGlobalEquations-1) NumEntries--;
     assert(B.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values+Valstart, Indices)==0);
    }
  
  // Finish up
  assert(B.TransformToLocal()==0);
  if (NumGlobalEquations<20) cout << "\n\nMatrix B = " << B << endl;


  Petra_Time & timer = *new Petra_Time(Comm);
  assert(B.Multiply(false, Y, BY)==0); // Compute BY = B*Y
  double elapsed_time = timer.ElapsedTime();
  double total_flops = B.Flops();
  double MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) cout << "\n\nTotal MFLOPs for B*Y = " << MFLOPs << endl<< endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////

  timer.ResetStartTime();
  assert(A.Multiply(true, Y, X)==0); // Compute X = A^T*Y
  elapsed_time = timer.ElapsedTime();
  total_flops = A.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) cout << "\n\nTotal MFLOPs for A^T*Y = " << MFLOPs << endl<< endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////

  // Iterate
  timer.ResetStartTime();
  assert(A.Multiply(false, X, AATY)==0); // Compute AATY = A*X
  elapsed_time = timer.ElapsedTime();
  total_flops = A.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;
  Petra_RDP_Vector & resid = *new Petra_RDP_Vector(YMap);
  resid.Update(1.0, BY, -1.0, AATY, 0.0);
  double residual;
  resid.Norm2(&residual);

  if (verbose) cout << "\n\nTotal MFLOPs for A*X = " << MFLOPs << endl<< endl;
  if (verbose) cout << "Residual = " << residual << endl<< endl;

  // Release all objects
  delete &resid;
  delete [] Values;
  delete [] Indices;
  delete [] MyGlobalElements;
  delete [] XGlobalElements;
  delete [] YGlobalElements;

  delete &timer;
  delete &X;
  delete &BY;
  delete &Y;
  delete &AATY;
  delete &XMap;
  delete &A;
  delete &RowMap;
			

  delete &Comm;
			
#ifdef PETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}

