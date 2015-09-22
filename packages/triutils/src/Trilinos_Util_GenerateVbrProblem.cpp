// @HEADER
// ***********************************************************************
// 
//                 TriUtils: Trilinos Utilities Package
//                 Copyright (2011) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
// @HEADER

#include "Trilinos_Util.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_VbrMatrix.h"

// CJ TODO FIXME: Trilinos_Util_GenerateVbrProblem available only if 32 bit GIDs available.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES

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

// npoints (In) - number of points in finite difference stencil
// xoff    (In) - stencil offsets in x direction (of length npoints)
// yoff    (In) - stencil offsets in y direction (of length npoints)
//   A standard 5-point finite difference stencil would be described as:
//     npoints = 5
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
// x      (Out) - Initial guess vector set to zero.
// b      (Out) - Generated RHS.  Values satisfy b = A*xexact
// xexact (Out) - Generated exact solution to Ax = b.

// Note: Caller of this function is responsible for deleting all output objects.

void Trilinos_Util_GenerateVbrProblem(int nx, int ny, int npoints, int * xoff, int * yoff,
																			int nsizes, int * sizes,
																			const Epetra_Comm  &comm, 
																			Epetra_BlockMap *& map, 
																			Epetra_VbrMatrix *& A, 
																			Epetra_Vector *& x, 
																			Epetra_Vector *& b,
																			Epetra_Vector *&xexact) {
	
	Epetra_MultiVector * x1, * b1, * xexact1;
	
	Trilinos_Util_GenerateVbrProblem(nx, ny, npoints, xoff, yoff, nsizes, sizes,
																	 1, comm, map, A, x1, b1, xexact1);

	x = dynamic_cast<Epetra_Vector *>(x1);
	b = dynamic_cast<Epetra_Vector *>(b1);
	xexact = dynamic_cast<Epetra_Vector *>(xexact1);

	return;
}

void Trilinos_Util_GenerateVbrProblem(int nx, int ny, int npoints, int * xoff, int * yoff, 
																			int nsizes, int * sizes, int nrhs,
																			const Epetra_Comm  &comm, 
																			Epetra_BlockMap *& map, 
																			Epetra_VbrMatrix *& A, 
																			Epetra_MultiVector *& x, 
																			Epetra_MultiVector *& b,
																			Epetra_MultiVector *&xexact) {

	int i, j;

	// Number of global equations is nx*ny.  These will be distributed in a linear fashion
	int numGlobalEquations = nx*ny;
  Epetra_Map ptMap(numGlobalEquations, 0, comm); // Create map with equal distribution of equations.

	int numMyElements = ptMap.NumMyElements();

	Epetra_IntVector elementSizes(ptMap); // This vector will have the list of element sizes
	for (i=0; i<numMyElements; i++) 
		elementSizes[i] = sizes[ptMap.GID64(i)%nsizes]; // cycle through sizes array

	map = new Epetra_BlockMap(-1, numMyElements, ptMap.MyGlobalElements(), elementSizes.Values(),
														ptMap.IndexBase(), ptMap.Comm());

  
  A = new Epetra_VbrMatrix(Copy, *map, 0); // Construct matrix

	int * indices = new int[npoints];
//	double * values = new double[npoints];

//	double dnpoints = (double) npoints;

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

#endif
