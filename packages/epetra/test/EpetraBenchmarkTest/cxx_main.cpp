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

#include <iostream>
#include <cstring>
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Time.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Version.h"

// prototypes

void GenerateCrsProblem(int numNodesX, int numNodesY, int numProcsX, int numProcsY, int numPoints, 
			int * xoff, int * yoff,
			const Epetra_Comm  &comm, 
			Epetra_Map *& map, 
			Epetra_CrsMatrix *& A, 
			Epetra_Vector *& b, 
			Epetra_Vector *& bt,
			Epetra_Vector *&xexact, bool StaticProfile, bool MakeLocalOnly);

void GenerateCrsProblem(int numNodesX, int numNodesY, int numProcsX, int numProcsY, int numPoints, 
			int * xoff, int * yoff, int nrhs,
			const Epetra_Comm  &comm, 
			Epetra_Map *& map, 
			Epetra_CrsMatrix *& A, 
			Epetra_MultiVector *& b, 
			Epetra_MultiVector *& bt,
			Epetra_MultiVector *&xexact, bool StaticProfile, bool MakeLocalOnly);
 

void GenerateMyGlobalElements(int numNodesX, int numNodesY, int numProcsX, int numProcs,
			      int myPID, int * & myGlobalElements);

void runMatrixTests(Epetra_CrsMatrix * A,  Epetra_MultiVector * b, Epetra_MultiVector * bt,
		    Epetra_MultiVector * xexact, bool StaticProfile);

int main(int argc, char *argv[])
{

  const string EpetraBenchmarkTest = "Epetra Benchmark Test Version 0.2 08/30/2007"; // Change this number and date when changes are made to functionality
  int ierr = 0;
  double elapsed_time;
  double total_flops;
  double MFLOPs;
  double global_dimension;
  double global_nonzero_count;
    

#ifdef EPETRA_MPI

  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm comm;
#endif

  // Check if we should print only selected timing info to standard out
  bool printflops = true, printmflops = true, printtime = true;
  if (argc>5) if (argv[5][0]=='-' && argv[5][1]=='f') {printmflops = false; printtime = false;}
  bool mflopsonly = false;
  if (argc>5) if (argv[5][0]=='-' && argv[5][1]=='m') {printflops = false; printtime = false;}
  bool timeonly = false;
  if (argc>5) if (argv[5][0]=='-' && argv[5][1]=='t') {printflops = false; printmflops = false;}

  // Check if we should print header to standard out
  bool makeheader = false;
  if (argc>6) if (argv[6][0]=='-' && argv[6][1]=='h') makeheader = true;

  if(argc < 5) {
    cerr << "Usage: " << argv[0]
         << " NumNodesX NumNodesY NumProcX NumProcY [-a|-f|-m|-t [-h]]" << endl
         << "where:" << endl
         << "NumNodesX         - Number of mesh nodes in X direction per processor" << endl
         << "NumNodesY         - Number of mesh nodes in Y direction per processor" << endl
         << "NumProcX          - Number of processors to use in X direction" << endl
         << "NumProcY          - Number of processors to use in Y direction" << endl
	 << "-a|-f|-m|-t       - Type of information to print: a=all, f=FLOPS, m=MFLOP/s, t=time (sec).  Default is -a."
         << "-h                - (Optional) Printer output table header if -h present (typically done on first call)." << endl
         << " NOTES: NumProcX*NumProcY must equal the number of processors used to run the problem." << endl << endl
	 << " Serial example:" << endl << endl
         << argv[0] << " 200 300 1 1 -m" << endl
	 << " Run this program on 1 processor using a 200 by 300 grid, printing only MFLOP/s information without a header."<< endl <<endl
	 << " MPI example:" << endl << endl
         << "mpirun -np 32 " << argv[0] << " 250 200 4 8 -a -h" << endl
	 << " Run this program on 32 processors putting a 250 by 200 subgrid on each processor using 4 processors "<< endl
	 << " in the X direction and 8 in the Y direction.  Total grid size is 1000 points in X and 1600 in Y for a total of 1.6M equations."<< endl
	 << " Print all information. Print header." << endl
         << endl;
    return(1);

  }
    //char tmp;
    //if (comm.MyPID()==0) cout << "Press any key to continue..."<< endl;
    //if (comm.MyPID()==0) cin >> tmp;
    //comm.Barrier();

  if (makeheader && comm.MyPID()==0)
    cout << EpetraBenchmarkTest << endl
	 << "Using " << Epetra_Version() << endl << endl;
  if (makeheader) cout << comm <<endl;


  // Redefine makeheader to only print on PE 0

  if (makeheader && comm.MyPID()!=0) makeheader = false;

  int numNodesX = atoi(argv[1]);
  int numNodesY = atoi(argv[2]);
  int numProcsX = atoi(argv[3]);
  int numProcsY = atoi(argv[4]);
  int numPoints = 25;

  if (makeheader) {
    cout << " Number of local nodes in X direction  = " << numNodesX << endl
	 << " Number of local nodes in Y direction  = " << numNodesY << endl
	 << " Number of global nodes in X direction = " << numNodesX*numProcsX << endl
	 << " Number of global nodes in Y direction = " << numNodesY*numProcsY << endl
	 << " Number of local nonzero entries       = " << numNodesX*numNodesY*numPoints << endl
	 << " Number of global nonzero entries      = " << numNodesX*numNodesY*numPoints*numProcsX*numProcsY << endl
	 << " Number of Processors in X direction   = " << numProcsX << endl
	 << " Number of Processors in Y direction   = " << numProcsY << endl
	 << " Number of Points in stencil           = " << numPoints << endl << endl;
    cout << " Timing the following:" <<endl
	 << " SpMV - Sparse matrix vector product using Epetra_CrsMatrix class" << endl
	 << " SpMM2- Sparse matrix times 2-column multivector using Epetra_CrsMatrix class" << endl
	 << " SpMM4- Sparse matrix times 4-column multivector using Epetra_CrsMatrix class" << endl
	 << " SpMM8- Sparse matrix times 8-column multivector using Epetra_CrsMatrix class" << endl
	 << " 2-norm of an Epetra_MultiVector" << endl
	 << " Dot-product of 2 Epetra_MultiVectors" << endl
	 << " AXPY of 2 Epetra_MultiVectors" << endl << endl;
  }

  if (numProcsX*numProcsY!=comm.NumProc()) {
    cerr << "Number of processors = " << comm.NumProc() << endl
	 << " is not the product of " << numProcsX << " and " << numProcsY << endl << endl;
    return(1);
  }

  if (numNodesX*numNodesY<=0) {
    cerr << "Product of number of nodes is <= zero" << endl << endl;
    return(1);
  }

  Epetra_IntSerialDenseVector Xoff, XLoff, XUoff;
  Epetra_IntSerialDenseVector Yoff, YLoff, YUoff;
  // Generate a 25-point 2D Finite Difference matrix
  Xoff.Size(25);
  Yoff.Size(25);
  int xi = 0, yi = 0;
  int xo = -2, yo = -2;
  Xoff[xi++] = xo++;  Xoff[xi++] = xo++; Xoff[xi++] = xo++; Xoff[xi++] = xo++; Xoff[xi++] = xo++;
  Yoff[yi++] = yo  ;  Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; 
  xo = -2, yo++;
  Xoff[xi++] = xo++;  Xoff[xi++] = xo++; Xoff[xi++] = xo++; Xoff[xi++] = xo++; Xoff[xi++] = xo++;
  Yoff[yi++] = yo  ;  Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; 
  xo = -2, yo++;
  Xoff[xi++] = xo++;  Xoff[xi++] = xo++; Xoff[xi++] = xo++; Xoff[xi++] = xo++; Xoff[xi++] = xo++;
  Yoff[yi++] = yo  ;  Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; 
  xo = -2, yo++;
  Xoff[xi++] = xo++;  Xoff[xi++] = xo++; Xoff[xi++] = xo++; Xoff[xi++] = xo++; Xoff[xi++] = xo++;
  Yoff[yi++] = yo  ;  Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; 
  xo = -2, yo++;
  Xoff[xi++] = xo++;  Xoff[xi++] = xo++; Xoff[xi++] = xo++; Xoff[xi++] = xo++; Xoff[xi++] = xo++;
  Yoff[yi++] = yo  ;  Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; Yoff[yi++] = yo  ; 
  
  

  Epetra_Map * map;
  Epetra_CrsMatrix * A;
  Epetra_MultiVector * b;
  Epetra_MultiVector * bt;
  Epetra_MultiVector * xexact;
  Epetra_SerialDenseVector resvec(0);
  global_dimension = (double) (numNodesX*numNodesY); // local grid dimension
  global_dimension *= (double) (numProcsX*numProcsY); // times number of processors
  global_nonzero_count = global_dimension * (double) numPoints; // number of nonzeros (w/o accounting for boundaries)

  //Timings
  Epetra_Time timer(comm);
  const int maxtimings = 7;
  const int maxtrials = 10;
  double results[maxtimings][3]; // timing results stored here

  int nrhs = 1;
  int ntimings = 0;
  for (int k=0; k<4; k++) {
    if (k>0) nrhs = nrhs*2;
      
    GenerateCrsProblem(numNodesX, numNodesY, numProcsX, numProcsY, numPoints,
		       Xoff.Values(), Yoff.Values(), nrhs, comm,
		       map, A, b, bt, xexact, true, false);

      
    Epetra_MultiVector z(*b);
    Epetra_MultiVector r(*b);
    Epetra_SerialDenseVector resvec(b->NumVectors());
    
    //Timings
    Epetra_Time timer(A->Comm());
    A->OptimizeStorage();
    
    timer.ResetStartTime();
    
    //maxtrials matvecs
    for( int i = 0; i < maxtrials; ++i )
      A->Multiply(false, *xexact, z); // Compute z = A*xexact or z = A'*xexact
    
    double elapsed_time = timer.ElapsedTime();
    double total_flops = 2.0*global_nonzero_count *((double) maxtrials);
    
    // Compute residual
    r.Update(-1.0, z, 1.0, *b, 0.0); // r = b - z
    
    r.Norm2(resvec.Values());
    double diff = resvec.NormInf();
    if (diff>1.0e-8 && comm.MyPID()==0) cerr << "Warning: Residual vector unusually large = " << diff << endl;

    double MFLOPs = total_flops/elapsed_time/1000000.0;

    results[ntimings][0] = total_flops;
    results[ntimings][1] = elapsed_time;
    results[ntimings++][2] = MFLOPs;

    delete A;
    delete b;
    delete bt; 
    delete xexact;
  } // end of k loop


  // *************** Vector timings *********************

  if (ntimings+3>maxtimings) cerr << "Variable maxtimings = " << maxtimings << " must be at least = " << ntimings+3 << endl;

  Epetra_MultiVector q(*map, nrhs);
  Epetra_MultiVector z(q);
  Epetra_MultiVector r(q);
  
  delete map;

  resvec.Resize(nrhs);
  
  
  timer.ResetStartTime();
  
  //maxtrials norms
  for( int i = 0; i < maxtrials; ++i )
    q.Norm2( resvec.Values() );
  
  elapsed_time = timer.ElapsedTime();
  total_flops = 2.0*global_dimension *((double) maxtrials);
  MFLOPs = total_flops/elapsed_time/1000000.0;
  results[ntimings][0] = total_flops;
  results[ntimings][1] = elapsed_time;
  results[ntimings++][2] = MFLOPs;



  timer.ResetStartTime();
  
  //maxtrials dot's
  for( int i = 0; i < maxtrials; ++i )
    q.Dot(z, resvec.Values());
  
  elapsed_time = timer.ElapsedTime();
  total_flops = 2.0*global_dimension *((double) maxtrials);
  MFLOPs = total_flops/elapsed_time/1000000.0;
  results[ntimings][0] = total_flops;
  results[ntimings][1] = elapsed_time;
  results[ntimings++][2] = MFLOPs;
  
  timer.ResetStartTime();
  
  //maxtrials updates
  for( int i = 0; i < maxtrials; ++i )
    q.Update(1.0, z, 1.0, r, 0.0);
  
  elapsed_time = timer.ElapsedTime();
  total_flops = 2.0*global_dimension *((double) maxtrials);
  MFLOPs = total_flops/elapsed_time/1000000.0;
  results[ntimings][0] = total_flops;
  results[ntimings][1] = elapsed_time;
  results[ntimings++][2] = MFLOPs;

  if (makeheader) 
    cout << "Metric_\t\t_Procs_\t__SpMV__\t_SpMM2__\t_SpMM4__\t_SpMM8__\t__NORM___\t__DOT___\t__AXPY__" << endl;
  if (comm.MyPID()==0) {
    cout.setf(std::ios::scientific);
    cout.precision(2);
    for (int j=0; j<3; j++) {
      bool doloop = false;
      if (j==0 && printflops) {	cout << "FLOPS\t"; doloop = true;}
      else if (j==1 && printtime) {cout << "Time(s)\t"; doloop = true;}
      else if (j==2 && printmflops) {cout << "MFLOP/s\t"; doloop = true;}
      if (doloop) {
	cout << "\t" << comm.NumProc();
	for (int i=0; i<maxtimings; i++)
	  cout << "\t" << results[i][j];
	cout << endl;
      }
    }
  }

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

return ierr ;
}

// Constructs a 2D PDE finite difference matrix using the list of x and y offsets.
// 
// nx      (In) - number of grid points in x direction
// ny      (In) - number of grid points in y direction
//   The total number of equations will be nx*ny ordered such that the x direction changes
//   most rapidly: 
//      First equation is at point (0,0)
//      Second at                  (1,0)
//       ...
//      nx equation at             (nx-1,0)
//      nx+1st equation at         (0,1)

// numPoints (In) - number of points in finite difference stencil
// xoff    (In) - stencil offsets in x direction (of length numPoints)
// yoff    (In) - stencil offsets in y direction (of length numPoints)
//   A standard 5-point finite difference stencil would be described as:
//     numPoints = 5
//     xoff = [-1, 1, 0,  0, 0]
//     yoff = [ 0, 0, 0, -1, 1]

// nrhs - Number of rhs to generate. (First interface produces vectors, so nrhs is not needed

// comm    (In) - an Epetra_Comm object describing the parallel machine (numProcs and my proc ID)
// map    (Out) - Epetra_Map describing distribution of matrix and vectors/multivectors
// A      (Out) - Epetra_CrsMatrix constructed for nx by ny grid using prescribed stencil
//                Off-diagonal values are random between 0 and 1.  If diagonal is part of stencil,
//                diagonal will be slightly diag dominant.
// b      (Out) - Generated RHS.  Values satisfy b = A*xexact
// bt     (Out) - Generated RHS.  Values satisfy b = A'*xexact
// xexact (Out) - Generated exact solution to Ax = b and b' = A'xexact

// Note: Caller of this function is responsible for deleting all output objects.

void GenerateCrsProblem(int numNodesX, int numNodesY, int numProcsX, int numProcsY, int numPoints, 
			int * xoff, int * yoff,
			const Epetra_Comm  &comm, 
			Epetra_Map *& map, 
			Epetra_CrsMatrix *& A, 
			Epetra_Vector *& b, 
			Epetra_Vector *& bt,
			Epetra_Vector *&xexact, bool StaticProfile, bool MakeLocalOnly) {

  Epetra_MultiVector * b1, * bt1, * xexact1;
	
  GenerateCrsProblem(numNodesX, numNodesY, numProcsX, numProcsY, numPoints, 
		     xoff, yoff, 1, comm, 
		     map, A, b1, bt1, xexact1, StaticProfile, MakeLocalOnly);

  b = dynamic_cast<Epetra_Vector *>(b1);
  bt = dynamic_cast<Epetra_Vector *>(bt1);
  xexact = dynamic_cast<Epetra_Vector *>(xexact1);

  return;
}

void GenerateCrsProblem(int numNodesX, int numNodesY, int numProcsX, int numProcsY, int numPoints, 
			int * xoff, int * yoff, int nrhs,
			const Epetra_Comm  &comm,
			Epetra_Map *& map, 
			Epetra_CrsMatrix *& A, 
			Epetra_MultiVector *& b, 
			Epetra_MultiVector *& bt,
			Epetra_MultiVector *&xexact, bool StaticProfile, bool MakeLocalOnly) {
  
  Epetra_Time timer(comm);
  // Determine my global IDs
  int * myGlobalElements;
  GenerateMyGlobalElements(numNodesX, numNodesY, numProcsX, numProcsY, comm.MyPID(), myGlobalElements);

  int numMyEquations = numNodesX*numNodesY;
  
  map = new Epetra_Map(-1, numMyEquations, myGlobalElements, 0, comm); // Create map with 2D block partitioning.
  delete [] myGlobalElements;

  int numGlobalEquations = map->NumGlobalElements();

  int profile = 0; if (StaticProfile) profile = numPoints;

#ifdef EPETRA_HAVE_STATICPROFILE

  if (MakeLocalOnly) 
    A = new Epetra_CrsMatrix(Copy, *map, *map, profile, StaticProfile); // Construct matrix with rowmap=colmap
  else
    A = new Epetra_CrsMatrix(Copy, *map, profile, StaticProfile); // Construct matrix

#else

  if (MakeLocalOnly) 
    A = new Epetra_CrsMatrix(Copy, *map, *map, profile); // Construct matrix with rowmap=colmap
  else
    A = new Epetra_CrsMatrix(Copy, *map, profile); // Construct matrix

#endif

  int * indices = new int[numPoints];
  double * values = new double[numPoints];

  double dnumPoints = (double) numPoints;
  int nx = numNodesX*numProcsX;

  for (int i=0; i<numMyEquations; i++) {

    int rowID = map->GID(i);
    int numIndices = 0;

    for (int j=0; j<numPoints; j++) {
      int colID = rowID + xoff[j] + nx*yoff[j]; // Compute column ID based on stencil offsets
      if (colID>-1 && colID<numGlobalEquations) {
	indices[numIndices] = colID;
	double value = - ((double) rand())/ ((double) RAND_MAX);
	if (colID==rowID)
	  values[numIndices++] = dnumPoints - value; // Make diagonal dominant
	else
	  values[numIndices++] = value;
      }
    }
    //cout << "Building row " << rowID << endl;
    A->InsertGlobalValues(rowID, numIndices, values, indices);
  }

  delete [] indices;
  delete [] values;
  double insertTime = timer.ElapsedTime();
  timer.ResetStartTime();
  A->FillComplete();
  double fillCompleteTime = timer.ElapsedTime();

  if (nrhs<=1) {  
    b = new Epetra_Vector(*map);
    bt = new Epetra_Vector(*map);
    xexact = new Epetra_Vector(*map);
  }
  else {
    b = new Epetra_MultiVector(*map, nrhs);
    bt = new Epetra_MultiVector(*map, nrhs);
    xexact = new Epetra_MultiVector(*map, nrhs);
  }

  xexact->Random(); // Fill xexact with random values

  A->Multiply(false, *xexact, *b);
  A->Multiply(true, *xexact, *bt);

  return;
}




void GenerateMyGlobalElements(int numNodesX, int numNodesY, int numProcsX, int numProcs,
			      int myPID, int * & myGlobalElements) {

  myGlobalElements = new int[numNodesX*numNodesY];
  int myProcX = myPID%numProcsX;
  int myProcY = myPID/numProcsX;
  int curGID = myProcY*(numProcsX*numNodesX)*numNodesY+myProcX*numNodesX;
  for (int j=0; j<numNodesY; j++) {
    for (int i=0; i<numNodesX; i++) {
      myGlobalElements[j*numNodesX+i] = curGID+i;
    }
    curGID+=numNodesX*numProcsX;
  }
  //for (int i=0; i<numNodesX*numNodesY; i++) cout << "MYPID " << myPID <<" GID "<< myGlobalElements[i] << endl;
  
  return;
}

