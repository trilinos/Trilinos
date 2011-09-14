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


#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Time.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Flops.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "../epetra_test_err.h"
#include "Epetra_IntVector.h"
#include "Epetra_Version.h"
#include "Epetra_RowMatrixTransposer.h"
#include "Epetra_Time.h"

int checkResults(Epetra_RowMatrix * A,
                 Epetra_CrsMatrix * transA,
                 Epetra_MultiVector * xexact,
                 bool verbose);

void GenerateCrsProblem(int nx, int ny, int npoints,
                        int * xoff, int * yoff, int nrhs,
                        const Epetra_Comm  &comm,
                        Epetra_Map *& map,
                        Epetra_CrsMatrix *& A,
                        Epetra_MultiVector *& x,
                        Epetra_MultiVector *& b,
                        Epetra_MultiVector *&xexact);

void GenerateVbrProblem(int nx, int ny, int npoints, int * xoff, int * yoff,
                        int nsizes, int * sizes, int nrhs,
                        const Epetra_Comm  &comm,
                        Epetra_BlockMap *& map,
                        Epetra_VbrMatrix *& A,
                        Epetra_MultiVector *& x,
                        Epetra_MultiVector *& b,
                        Epetra_MultiVector *&xexact);

int main(int argc, char *argv[]) {

  int ierr = 0, i;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);

  Epetra_MpiComm Comm( MPI_COMM_WORLD );

#else

  Epetra_SerialComm Comm;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  int verbose_int = verbose ? 1 : 0;
  Comm.Broadcast(&verbose_int, 1, 0);
  verbose = verbose_int==1 ? true : false;


  //  char tmp;
  //  if (Comm.MyPID()==0) cout << "Press any key to continue..."<< endl;
  //  if (Comm.MyPID()==0) cin >> tmp;
  //  Comm.Barrier();

  Comm.SetTracebackMode(0); // This should shut down any error traceback reporting
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  if(verbose && MyPID==0)
    cout << Epetra_Version() << endl << endl;

  int nx = 128;
  int ny = NumProc*nx; // Scale y grid with number of processors

  // Create funky stencil to make sure the matrix is non-symmetric (transpose non-trivial):

  // (i-1,j-1) (i-1,j  )
  // (i  ,j-1) (i  ,j  ) (i  ,j+1)
  // (i+1,j-1) (i+1,j  )

  int npoints = 7;

  int xoff[] = {-1,  0,  1, -1,  0,  1,  0};
  int yoff[] = {-1, -1, -1,  0,  0,  0,  1};

  Epetra_Map * map;
  Epetra_CrsMatrix * A;
  Epetra_MultiVector * x, * b, * xexact;

  GenerateCrsProblem(nx, ny, npoints, xoff, yoff, 1, Comm, map, A, x, b, xexact);

  if (nx<8) {
    cout << *A << endl;
    cout << "X exact = " << endl << *xexact << endl;
    cout << "B       = " << endl << *b << endl;
  }
  // Construct transposer object

  Epetra_Time timer(Comm);

  double start = timer.ElapsedTime();
  Epetra_RowMatrixTransposer transposer(A);
  if (verbose) cout << "\nTime to construct transposer  = "
       << timer.ElapsedTime() - start << endl;

  bool MakeDataContiguous = true;
  Epetra_CrsMatrix * transA;

  start = timer.ElapsedTime();
  transposer.CreateTranspose(MakeDataContiguous, transA);
  if (verbose) cout << "\nTime to create transpose matrix  = "
       << timer.ElapsedTime() - start << endl;


  // Now test output of transposer by performing matvecs
  ierr += checkResults(A, transA, xexact, verbose);


  // Now change values in original matrix and test update facility of transposer
  // Add 2 to the diagonal of each row

  double Value = 2.0;
  for (i=0; i< A->NumMyRows(); i++)
    A->SumIntoMyValues(i, 1, &Value, &i);


  start = timer.ElapsedTime();
  transposer.UpdateTransposeValues(A);
  if (verbose) cout << "\nTime to update transpose matrix  = "
     << timer.ElapsedTime() - start << endl;

  ierr += checkResults(A, transA, xexact, verbose);

  delete A;
  delete b;
  delete x;
  delete xexact;
  delete map;


  if (verbose) cout << endl << "Checking transposer for VbrMatrix objects" << endl<< endl;

  int nsizes = 4;
  int sizes[] = {4, 6, 5, 3};

  Epetra_VbrMatrix * Avbr;
  Epetra_BlockMap * bmap;

  GenerateVbrProblem(nx, ny, npoints, xoff, yoff, nsizes, sizes,
                     1,  Comm, bmap, Avbr, x, b, xexact);

  if (nx<8) {
    cout << *Avbr << endl;
    cout << "X exact = " << endl << *xexact << endl;
    cout << "B       = " << endl << *b << endl;
  }

  start = timer.ElapsedTime();
  Epetra_RowMatrixTransposer transposer1(Avbr);
  if (verbose) cout << "\nTime to construct transposer  = "
       << timer.ElapsedTime() - start << endl;

  start = timer.ElapsedTime();
  transposer1.CreateTranspose(MakeDataContiguous, transA);
  if (verbose) cout << "\nTime to create transpose matrix  = "
       << timer.ElapsedTime() - start << endl;


  // Now test output of transposer by performing matvecs
  ierr += checkResults(Avbr, transA, xexact, verbose);

  // Now change values in original matrix and test update facility of transposer
  // Scale matrix on the left by rowsums

  Epetra_Vector invRowSums(Avbr->RowMap());

  Avbr->InvRowSums(invRowSums);
  Avbr->LeftScale(invRowSums);

  start = timer.ElapsedTime();
  transposer1.UpdateTransposeValues(Avbr);
  if (verbose) cout << "\nTime to update transpose matrix  = "
    << timer.ElapsedTime() - start << endl;

  ierr += checkResults(Avbr, transA, xexact, verbose);

  delete Avbr;
  delete b;
  delete x;
  delete xexact;
  delete bmap;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

/* end main */
  return ierr ;
}

int checkResults(Epetra_RowMatrix * A, Epetra_CrsMatrix * transA,
                 Epetra_MultiVector * xexact, bool verbose)
{
  int n = A->NumGlobalRows();

  if (n<100) cout << "A transpose = " << endl << *transA << endl;

  Epetra_MultiVector x1(View, A->OperatorRangeMap(), &((*xexact)[0]), 1);
  Epetra_MultiVector b1(A->OperatorDomainMap(), 1);

  A->SetUseTranspose(true);

  Epetra_Time timer(A->Comm());
  double start = timer.ElapsedTime();
  A->Apply(x1, b1);
  if (verbose) cout << "\nTime to compute b1: matvec with original matrix using transpose flag  = " << timer.ElapsedTime() - start << endl;

  if (n<100) cout << "b1 = " << endl << b1 << endl;
  Epetra_MultiVector x2(View, transA->OperatorDomainMap(), &((*xexact)[0]), 1);
  Epetra_MultiVector b2(transA->OperatorRangeMap(), 1);
  start = timer.ElapsedTime();
  transA->Multiply(false, x2, b2);
  if (verbose) cout << "\nTime to compute b2: matvec with transpose matrix                      = " << timer.ElapsedTime() - start << endl;

  if (n<100) cout << "b1 = " << endl << b1 << endl;

  double residual;
  Epetra_MultiVector resid(A->OperatorRangeMap(), 1);

  resid.Update(1.0, b1, -1.0, b2, 0.0);
  int ierr0 = resid.Norm2(&residual);
  assert(ierr0==0);
  if (verbose) cout << "Norm of b1 - b2 = " << residual << endl;

  int ierr = 0;

  if (residual > 1.0e-10) ierr++;

  if (ierr!=0 && verbose) {
    cerr << "Status: Test failed" << endl;
  }
  else if (verbose) cerr << "Status: Test passed" << endl;

  return(ierr);
}

void GenerateCrsProblem(int nx, int ny, int npoints,
                        int * xoff, int * yoff, int nrhs,
                        const Epetra_Comm  &comm,
                        Epetra_Map *& map,
                        Epetra_CrsMatrix *& A,
                        Epetra_MultiVector *& x,
                        Epetra_MultiVector *& b,
                        Epetra_MultiVector *&xexact)
{
  // Number of global equations is nx*ny.  These will be distributed in a linear fashion
  int numGlobalEquations = nx*ny;
  map = new Epetra_Map(numGlobalEquations, 0, comm); // Create map with equal distribution of equations.

  int numMyEquations = map->NumMyElements();

  A = new Epetra_CrsMatrix(Copy, *map, 0); // Construct matrix

  int * indices = new int[npoints];
  double * values = new double[npoints];

  double dnpoints = (double) npoints;

  for (int i=0; i<numMyEquations; i++) {

    int rowID = map->GID(i);
    int numIndices = 0;

    for (int j=0; j<npoints; j++) {
      int colID = rowID + xoff[j] + nx*yoff[j]; // Compute column ID based on stencil offsets
      if (colID>-1 && colID<numGlobalEquations) {
        indices[numIndices] = colID;
        double value = - ((double) rand())/ ((double) RAND_MAX);
        if (colID==rowID)
          values[numIndices++] = dnpoints - value; // Make diagonal dominant
        else
          values[numIndices++] = -value;
      }
    }

    A->InsertGlobalValues(rowID, numIndices, values, indices);
  }

  delete [] indices;
  delete [] values;

  A->FillComplete();

  if (nrhs<=1) {
    x = new Epetra_Vector(*map);
    b = new Epetra_Vector(*map);
    xexact = new Epetra_Vector(*map);
  }
  else {
    x = new Epetra_MultiVector(*map, nrhs);
    b = new Epetra_MultiVector(*map, nrhs);
    xexact = new Epetra_MultiVector(*map, nrhs);
  }

  xexact->Random(); // Fill xexact with random values

  A->Multiply(false, *xexact, *b);

  return;
}

void GenerateVbrProblem(int nx, int ny, int npoints, int * xoff, int * yoff,
                        int nsizes, int * sizes, int nrhs,
                        const Epetra_Comm  &comm,
                        Epetra_BlockMap *& map,
                        Epetra_VbrMatrix *& A,
                        Epetra_MultiVector *& x,
                        Epetra_MultiVector *& b,
                        Epetra_MultiVector *&xexact)
{
  int i, j;

  // Number of global equations is nx*ny.  These will be distributed in a linear fashion
  int numGlobalEquations = nx*ny;
  Epetra_Map ptMap(numGlobalEquations, 0, comm); // Create map with equal distribution of equations.

  int numMyElements = ptMap.NumMyElements();

  Epetra_IntVector elementSizes(ptMap); // This vector will have the list of element sizes
  for (i=0; i<numMyElements; i++)
    elementSizes[i] = sizes[ptMap.GID(i)%nsizes]; // cycle through sizes array

  map = new Epetra_BlockMap(-1, numMyElements, ptMap.MyGlobalElements(), elementSizes.Values(), ptMap.IndexBase(), ptMap.Comm());


  A = new Epetra_VbrMatrix(Copy, *map, 0); // Construct matrix

  int * indices = new int[npoints];

  // This section of code creates a vector of random values that will be used to create
  // light-weight dense matrices to pass into the VbrMatrix construction process.

  int maxElementSize = 0;
  for (i=0; i< nsizes; i++) maxElementSize = EPETRA_MAX(maxElementSize, sizes[i]);

  Epetra_LocalMap lmap(maxElementSize*maxElementSize, ptMap.IndexBase(), ptMap.Comm());
  Epetra_Vector randvec(lmap);
  randvec.Random();
  randvec.Scale(-1.0); // Make value negative


  for (i=0; i<numMyElements; i++) {
    int rowID = map->GID(i);
    int numIndices = 0;
    int rowDim = sizes[rowID%nsizes];
    for (j=0; j<npoints; j++) {
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

  A->FillComplete();

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
    x = new Epetra_Vector(*map);
    b = new Epetra_Vector(*map);
    xexact = new Epetra_Vector(*map);
  }
  else {
    x = new Epetra_MultiVector(*map, nrhs);
    b = new Epetra_MultiVector(*map, nrhs);
    xexact = new Epetra_MultiVector(*map, nrhs);
  }

  xexact->Random(); // Fill xexact with random values

  A->Multiply(false, *xexact, *b);

  return;
}

