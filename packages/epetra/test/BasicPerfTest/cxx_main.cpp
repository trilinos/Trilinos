//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
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
// ************************************************************************
//@HEADER

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
#include "Epetra_Flops.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "../epetra_test_err.h"
#include "Epetra_Version.h"
                                            
// prototypes

void GenerateCrsProblem(int numNodesX, int numNodesY, int numProcsX, int numProcsY, int numPoints, 
			int * xoff, int * yoff,
			const Epetra_Comm  &comm, 
			Epetra_Map *& map, 
			Epetra_CrsMatrix *& A, 
			Epetra_Vector *& b, 
			Epetra_Vector *& bt,
			Epetra_Vector *&xexact);

void GenerateCrsProblem(int numNodesX, int numNodesY, int numProcsX, int numProcsY, int numPoints, 
			int * xoff, int * yoff, int nrhs,
			const Epetra_Comm  &comm, 
			Epetra_Map *& map, 
			Epetra_CrsMatrix *& A, 
			Epetra_MultiVector *& b, 
			Epetra_MultiVector *& bt,
			Epetra_MultiVector *&xexact);
 
void GenerateVbrProblem(int numNodesX, int numNodesY, int numProcsX, int numProcsY, int numPoints, 
			int * xoff, int * yoff,
			int nsizes, int * sizes,
			const Epetra_Comm  &comm, 
			Epetra_BlockMap *& map, 
			Epetra_VbrMatrix *& A, 
			Epetra_Vector *& b, 
			Epetra_Vector *& bt,
			Epetra_Vector *&xexact);

void GenerateVbrProblem(int numNodesX, int numNodesY, int numProcsX, int numProcsY, int numPoints, 
			int * xoff, int * yoff, 
			int nsizes, int * sizes, int nrhs,
			const Epetra_Comm  &comm, 
			Epetra_BlockMap *& map, 
			Epetra_VbrMatrix *& A, 
			Epetra_MultiVector *& b, 
			Epetra_MultiVector *& bt,
			Epetra_MultiVector *&xexact);

void GenerateMyGlobalElements(int numNodesX, int numNodesY, int numProcsX, int numProcs,
			      int myPID, int * & myGlobalElements);

int main(int argc, char *argv[])
{
  int ierr = 0, j;
  double elapsed_time;
  double total_flops;
  double MFLOPs;
    

#ifdef EPETRA_MPI

  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm comm;
#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>6) if (argv[6][0]=='-' && argv[6][1]=='v') verbose = true;

  if(argc < 6) {
    cerr << "Usage: " << argv[0]
         << " NumNodesX NumNodesY NumProcX NumProcY NumPoints [-v]" << endl
         << "where:" << endl
         << "NumNodesX         - Number of mesh nodes in X direction per processor" << endl
         << "NumNodesY         - Number of mesh nodes in Y direction per processor" << endl
         << "NumProcX          - Number of processors to use in X direction" << endl
         << "NumProcY          - Number of processors to use in Y direction" << endl
         << "NumPoints         - Number of points to use in stencil (5 or 9 only)" << endl
         << "-v                - (Optional) Run in verbose mode if present" << endl
         << " NOTES: NumProcX*NumProcY must equal the number of processors used to run the problem. Example:" << endl
         << "mpirun -np 32 " << argv[0] << " 10 12 4 8 -v" << endl
	 << " Run this program on 32 processors putting a 10 X 12 subgrid on each processor using 4 processors "<< endl
	 << " in the X direction and 8 in the Y direction.  Total grid size is 40 points in X and 96 in Y."<< endl
         << endl;
    return(1);

  }
  //  char tmp;
  //  if (rank==0) cout << "Press any key to continue..."<< endl;
  //  if (rank==0) cin >> tmp;
  //  comm.Barrier();

  comm.SetTracebackMode(0); // This should shut down any error traceback reporting
  if (verbose && comm.MyPID()==0)
    cout << Epetra_Version() << endl << endl;

  if (verbose) cout << comm <<endl;

  //  bool verbose1 = verbose;

  // Redefine verbose to only print on PE 0

  if (verbose && comm.MyPID()!=0) verbose = false;

  int numNodesX = atoi(argv[1]);
  int numNodesY = atoi(argv[2]);
  int numProcsX = atoi(argv[3]);
  int numProcsY = atoi(argv[4]);
  int numPoints = atoi(argv[5]);

  if (verbose) {
    cout << " Number of local nodes in X direction  = " << numNodesX << endl
	 << " Number of local nodes in Y direction  = " << numNodesY << endl
	 << " Number of global nodes in X direction = " << numNodesX*numProcsX << endl
	 << " Number of global nodes in Y direction = " << numNodesY*numProcsY << endl
	 << " Number of local nonzero entries       = " << numNodesX*numNodesY*numPoints << endl
	 << " Number of global nonzero entries      = " << numNodesX*numNodesY*numPoints*numProcsX*numProcsY << endl
	 << " Number of Processors in X direction   = " << numProcsX << endl
	 << " Number of Processors in Y direction   = " << numProcsY << endl
	 << " Number of Points in stencil           = " << numPoints << endl << endl;
  }

  if (numProcsX*numProcsY!=comm.NumProc()) {
    cerr << "Number of processors = " << comm.NumProc() << endl
	 << " is not the product of " << numProcsX << " and " << numProcsY << endl << endl;
    return(1);
  }

  if (numPoints!=5 && numPoints!=9) {
    cerr << "Number of points specified = " << numPoints << endl
	 << " is not 5 or 9" << endl << endl;
    return(1);
  }

  if (numNodesX*numNodesY<=0) {
    cerr << "Product of number of nodes is <= zero" << endl << endl;
    return(1);
  }

  Epetra_IntSerialDenseVector Xoff;
  Epetra_IntSerialDenseVector Yoff;
  if (numPoints==5) {

     // Generate a 5-point 2D Finite Difference matrix
    Xoff.Size(5);
    Yoff.Size(5);
    Xoff[0] = -1; Xoff[1] = 1; Xoff[2] = 0; Xoff[3] = 0;  Xoff[4] = 0; 
    Yoff[0] = 0;  Yoff[1] = 0; Yoff[2] = 0; Yoff[3] = -1; Yoff[4] = 1; 
  }
  else {
    // Generate a 9-point 2D Finite Difference matrix
    Xoff.Size(9);
    Yoff.Size(9);
    Xoff[0] = -1;  Xoff[1] =  0; Xoff[2] =  1; 
    Yoff[0] = -1;  Yoff[1] = -1; Yoff[2] = -1; 
    Xoff[3] = -1;  Xoff[4] =  0; Xoff[5] =  1; 
    Yoff[3] =  0;  Yoff[4] =  0; Yoff[5] =  0; 
    Xoff[6] = -1;  Xoff[7] =  0; Xoff[8] =  1; 
    Yoff[6] =  1;  Yoff[7] =  1; Yoff[8] =  1; 
  }

  Epetra_Map * map;
  Epetra_CrsMatrix * A;
  Epetra_Vector * b;
  Epetra_Vector * bt;
  Epetra_Vector * xexact;

  GenerateCrsProblem(numNodesX, numNodesY, numProcsX, numProcsY, numPoints,
 Xoff.Values(), Yoff.Values(), comm, 
		     map, A, b, bt, xexact);

  Epetra_Vector q(b->Map());
  Epetra_Vector z(b->Map());
  Epetra_Vector r(b->Map());

  //Timings
  Epetra_Flops flopcounter;
  A->SetFlopCounter(flopcounter);
  Epetra_Time timer(comm);
    
  for (j=0; j<2; j++) { // j = 0 is notrans, j = 1 is trans
      
    flopcounter.ResetFlops();
    timer.ResetStartTime();

    bool TransA = (j==1);
    //10 matvecs
    for( int i = 0; i < 10; ++i )
      A->Multiply(TransA, *xexact, z); // Compute z = A*xexact or z = A'*xexact
      
    elapsed_time = timer.ElapsedTime();
    total_flops = A->Flops();
    MFLOPs = total_flops/elapsed_time/1000000.0;
    if (verbose) cout << "\n\nTotal MFLOPs for 10 MatVec's (Trans = " << TransA
		      << ")     = " << MFLOPs << endl<< endl;
      

    // Compute residual
    if (TransA)
      r.Update(-1.0, z, 1.0, *bt, 0.0); // r = bt - z
    else
      r.Update(-1.0, z, 1.0, *b, 0.0); // r = b - z

    double rnorm;
    r.Norm2(&rnorm);
    if (verbose) cout << "Norm of difference between computed and exact RHS = " << rnorm << endl;
  }

  q.SetFlopCounter(*A);
  z.SetFlopCounter(*A);
    
  flopcounter.ResetFlops();
  timer.ResetStartTime();
  //10 norms
  double n_out;
  for( int i = 0; i < 10; ++i )
    q.Norm2( &n_out );

  elapsed_time = timer.ElapsedTime();
  total_flops = q.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;
  if (verbose) cout << "\n\nTotal MFLOPs for 10 Norm2's= " << MFLOPs << endl<< endl;
  r.SetFlopCounter(*A);
    
  flopcounter.ResetFlops();
  timer.ResetStartTime();

  //10 dot's
  for( int i = 0; i < 10; ++i )
    q.Dot(z, &n_out);
    
  elapsed_time = timer.ElapsedTime();
  total_flops = q.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;
  if (verbose) cout << "\n\nTotal MFLOPs for 10 Dot's= " << MFLOPs << endl<< endl;
    
  delete map;
  delete A;
  delete b;
  delete bt; 
  delete xexact;
		
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
			Epetra_Vector *&xexact) {

  Epetra_MultiVector * b1, * bt1, * xexact1;
	
  GenerateCrsProblem(numNodesX, numNodesY, numProcsX, numProcsY, numPoints, 
		     xoff, yoff, 1, comm, 
		     map, A, b1, bt1, xexact1);

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
			Epetra_MultiVector *&xexact) {
  
  Epetra_Time timer(comm);
  // Determine my global IDs
  int * myGlobalElements;
  GenerateMyGlobalElements(numNodesX, numNodesY, numProcsX, numProcsY, comm.MyPID(), myGlobalElements);

  int numMyEquations = numNodesX*numNodesY;
  
  map = new Epetra_Map(-1, numMyEquations, myGlobalElements, 0, comm); // Create map with 2D block partitioning.
  delete [] myGlobalElements;

  int numGlobalEquations = map->NumGlobalElements();
  
  A = new Epetra_CrsMatrix(Copy, *map, 0); // Construct matrix

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
	  values[numIndices++] = -value;
      }
    }
    //cout << "Building row " << rowID << endl;
    A->InsertGlobalValues(rowID, numIndices, values, indices);
  }

  delete [] indices;
  delete [] values;
  double insertTime = timer.ElapsedTime();
  timer.ResetStartTime();
  A->TransformToLocal();
  double fillCompleteTime = timer.ElapsedTime();

  if (comm.MyPID()==0)
    cout << "Time to insert matrix values = " << insertTime << endl
	 << "Time to complete fill        = " << fillCompleteTime << endl;

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

// nsizes  (In) - Length of element size list used to create variable block map and matrix
// sizes   (In) - integer list of element sizes of length nsizes
//    The map associated with this VbrMatrix will be created by cycling through the sizes list.
//    For example, if nsize = 3 and sizes = [ 2, 4, 3], the block map will have elementsizes
//    of 2, 4, 3, 2, 4, 3, ...

// nrhs - Number of rhs to generate. (First interface produces vectors, so nrhs is not needed

// comm    (In) - an Epetra_Comm object describing the parallel machine (numProcs and my proc ID)
// map    (Out) - Epetra_Map describing distribution of matrix and vectors/multivectors
// A      (Out) - Epetra_VbrMatrix constructed for nx by ny grid using prescribed stencil
//                Off-diagonal values are random between 0 and 1.  If diagonal is part of stencil,
//                diagonal will be slightly diag dominant.
// b      (Out) - Generated RHS.  Values satisfy b = A*xexact
// bt     (Out) - Generated RHS.  Values satisfy b = A'*xexact
// xexact (Out) - Generated exact solution to Ax = b and b' = A'xexact

// Note: Caller of this function is responsible for deleting all output objects.

void GenerateVbrProblem(int numNodesX, int numNodesY, int numProcsX, int numProcsY, int numPoints, 
			int * xoff, int * yoff,
			int nsizes, int * sizes,
			const Epetra_Comm  &comm, 
			Epetra_BlockMap *& map, 
			Epetra_VbrMatrix *& A, 
			Epetra_Vector *& b, 
			Epetra_Vector *& bt,
			Epetra_Vector *&xexact) {
	
  Epetra_MultiVector * b1, * bt1, * xexact1;
	
  GenerateVbrProblem(numNodesX, numNodesY, numProcsX, numProcsY, numPoints,
		     xoff, yoff, nsizes, sizes,
		     1, comm, map, A, b1, bt1, xexact1);

  b = dynamic_cast<Epetra_Vector *>(b1);
  bt = dynamic_cast<Epetra_Vector *>(bt1);
  xexact = dynamic_cast<Epetra_Vector *>(xexact1);

  return;
}

void GenerateVbrProblem(int numNodesX, int numNodesY, int numProcsX, int numProcsY, int numPoints, 
			int * xoff, int * yoff, 
			int nsizes, int * sizes, int nrhs,
			const Epetra_Comm  &comm, 
			Epetra_BlockMap *& map, 
			Epetra_VbrMatrix *& A, 
			Epetra_MultiVector *& b, 
			Epetra_MultiVector *& bt,
			Epetra_MultiVector *&xexact) {

  int i, j;

  // Determine my global IDs
  int * myGlobalElements;
  GenerateMyGlobalElements(numNodesX, numNodesY, numProcsX, numProcsY, comm.MyPID(), myGlobalElements);

  int numMyElements = numNodesX*numNodesY;
  
  Epetra_Map ptMap(-1, numMyElements, myGlobalElements, 0, comm); // Create map with 2D block partitioning.
  delete [] myGlobalElements;

  int numGlobalEquations = ptMap.NumGlobalElements();

  Epetra_IntVector elementSizes(ptMap); // This vector will have the list of element sizes
  for (i=0; i<numMyElements; i++) 
    elementSizes[i] = sizes[ptMap.GID(i)%nsizes]; // cycle through sizes array

  map = new Epetra_BlockMap(-1, numMyElements, ptMap.MyGlobalElements(), elementSizes.Values(),
			    ptMap.IndexBase(), ptMap.Comm());

  
  A = new Epetra_VbrMatrix(Copy, *map, 0); // Construct matrix

  int * indices = new int[numPoints];
  double * values = new double[numPoints];

  // This section of code creates a vector of random values that will be used to create
  // light-weight dense matrices to pass into the VbrMatrix construction process.

  int maxElementSize = 0;
  for (i=0; i< nsizes; i++) maxElementSize = EPETRA_MAX(maxElementSize, sizes[i]);

  Epetra_LocalMap lmap(maxElementSize*maxElementSize, ptMap.IndexBase(), ptMap.Comm());
  Epetra_Vector randvec(lmap);
  randvec.Random();
  randvec.Scale(-1.0); // Make value negative
  int nx = numNodesX*numProcsX;


  for (i=0; i<numMyElements; i++) {
    int rowID = map->GID(i);
    int numIndices = 0;
    int rowDim = sizes[rowID%nsizes];
    for (j=0; j<numPoints; j++) {
      int colID = rowID + xoff[j] + nx*yoff[j]; // Compute column ID based on stencil offsets
      if (colID>-1 && colID<numGlobalEquations)
	indices[numIndices++] = colID;
    }
			
    A->BeginInsertGlobalValues(rowID, numIndices, indices);
		
    for (j=0; j < numIndices; j++) {
      int colDim = sizes[indices[j]%nsizes];
      A->SubmitBlockEntry(&(randvec[0]), rowDim, rowDim, colDim);
    }
    A->EndSubmitEntries();
  }

  delete [] indices;

  A->TransformToLocal();

  // Compute the InvRowSums of the matrix rows
  Epetra_Vector invRowSums(A->RowMap());
  Epetra_Vector rowSums(A->RowMap());
  A->InvRowSums(invRowSums);
  rowSums.Reciprocal(invRowSums);

  // Jam the row sum values into the diagonal of the Vbr matrix (to make it diag dominant)
  int numBlockDiagonalEntries;
  int * rowColDims;
  int * diagoffsets = map->FirstPointInElementList();
  A->BeginExtractBlockDiagonalView(numBlockDiagonalEntries, rowColDims);
  for (i=0; i< numBlockDiagonalEntries; i++) {
    double * diagVals;
    int diagLDA;
    A->ExtractBlockDiagonalEntryView(diagVals, diagLDA);
    int rowDim = map->ElementSize(i);
    for (j=0; j<rowDim; j++) diagVals[j+j*diagLDA] = rowSums[diagoffsets[i]+j];
  }

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

