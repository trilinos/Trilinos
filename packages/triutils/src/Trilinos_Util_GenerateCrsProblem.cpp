// @HEADER
// ***********************************************************************
// 
//                 TriUtils: Trilinos Utilities Package
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
// @HEADER

#include "Trilinos_Util.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"

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

// nrhs - Number of rhs to generate. (First interface produces vectors, so nrhs is not needed

// comm    (In) - an Epetra_Comm object describing the parallel machine (numProcs and my proc ID)
// map    (Out) - Epetra_Map describing distribution of matrix and vectors/multivectors
// A      (Out) - Epetra_CrsMatrix constructed for nx by ny grid using prescribed stencil
//                Off-diagonal values are random between 0 and 1.  If diagonal is part of stencil,
//                diagonal will be slightly diag dominant.
// x      (Out) - Initial guess vector set to zero.
// b      (Out) - Generated RHS.  Values satisfy b = A*xexact
// xexact (Out) - Generated exact solution to Ax = b.

// Note: Caller of this function is responsible for deleting all output objects.

void Trilinos_Util_GenerateCrsProblem(int nx, int ny, int npoints, int * xoff, int * yoff,
				      const Epetra_Comm  &comm, 
				      Epetra_Map *& map, 
				      Epetra_CrsMatrix *& A, 
				      Epetra_Vector *& x, 
				      Epetra_Vector *& b,
				      Epetra_Vector *&xexact, int indexBase) {

  Epetra_MultiVector * x1, * b1, * xexact1;
	
  Trilinos_Util_GenerateCrsProblem(nx, ny, npoints, xoff, yoff, 1, comm, map, A, x1, b1, xexact1, indexBase);

  x = dynamic_cast<Epetra_Vector *>(x1);
  b = dynamic_cast<Epetra_Vector *>(b1);
  xexact = dynamic_cast<Epetra_Vector *>(xexact1);

  return;
}

void Trilinos_Util_GenerateCrsProblem(int nx, int ny, int npoints, int * xoff, int * yoff, int nrhs,
				      const Epetra_Comm  &comm, 
				      Epetra_Map *& map, 
				      Epetra_CrsMatrix *& A, 
				      Epetra_MultiVector *& x, 
				      Epetra_MultiVector *& b,
				      Epetra_MultiVector *&xexact, int indexBase) {

  // Number of global equations is nx*ny.  These will be distributed in a linear fashion
  int numGlobalEquations = nx*ny;
  map = new Epetra_Map(numGlobalEquations, indexBase, comm); // Create map with equal distribution of equations.
  
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
      if (colID>indexBase-1 && colID<indexBase+numGlobalEquations) {
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
