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

#include <iostream>
#include <string>
#include "Kokkos_DenseVector.hpp"
#include "Kokkos_DenseMultiVector.hpp"
#include "Kokkos_HbMatrix.hpp"

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

// A      (Out) - Kokkos::HbMatrix constructed for nx by ny grid using prescribed stencil
//                Off-diagonal values are random between 0 and 1.  If diagonal is part of stencil,
//                diagonal will be slightly diag dominant.
// x      (Out) - Initial guess vector set to zero.
// b      (Out) - Generated RHS.  Values satisfy b = A*xexact
// xexact (Out) - Generated exact solution to Ax = b.

// Note: Caller of this function is responsible for deleting all output objects.

template<typename OrdinalType, typename ScalarType>
void Trilinos_Util_GenerateCrsProblem(bool generateClassicHbMatrix, bool isRowOriented,
				      OrdinalType nx, OrdinalType ny, OrdinalType npoints, 
				      OrdinalType * xoff, OrdinalType * yoff,
				      Kokkos::HbMatrix<OrdinalType, ScalarType> *& A, 
				      Kokkos::Vector<OrdinalType, ScalarType> *& x, 
				      Kokkos::Vector<OrdinalType, ScalarType> *& b,
				      Kokkos::Vector<OrdinalType, ScalarType> *&xexact) {

  Kokkos::DenseMultiVector<OrdinalType, ScalarType> * x1, * b1, * xexact1;
	
  Trilinos_Util_GenerateCrsProblem(generateClassicHbmatrix, isRowOriented, nx, ny, npoints, xoff, yoff, 1, A, x1, b1, xexact1);

  x = new Kokkos::Vector<OrdinalType, ScalarType>();
  b = new Kokkos::Vector<OrdinalType, ScalarType>();
  xexact = new Kokkos::Vector<OrdinalType, ScalarType>();

  x.initializeValues(x1.getValues(0), x1.getValues(0), x1.getColInc());
  b.initializeValues(b1.getValues(0), b1.getValues(0), b1.getColInc());
  xexact.initializeValues(xexact1.getNumRows(), xexact1.getValues(0), xexact1.getColInc());

  return;
}

template<typename OrdinalType, typename ScalarType>
void Trilinos_Util_GenerateCrsProblem(bool generateClassicHbMatrix, bool isRowOriented,
				      OrdinalType nx, OrdinalType ny, OrdinalType npoints, 
				      OrdinalType * xoff, OrdinalType * yoff, OrdinalType nrhs,
				      Kokkos::HbMatrix<OrdinalType, ScalarType> *& A, 
				      Kokkos::MultiVector<OrdinalType, ScalarType> *& x, 
				      Kokkos::MultiVector<OrdinalType, ScalarType> *& b,
				      Kokkos::MultiVector<OrdinalType, ScalarType> *&xexact) {

  // Number of equations is nx*ny.
  OrdinalType numEquations = nx*ny;
  OrdinalType i, j, k;

  

    ScalarType ** xv = new ScalarType*[nrhs];
    ScalarType ** bv = new ScalarType*[nrhs];
    ScalarType ** xexactv = new ScalarType*[nrhs];

    for (i=0; i<nrhs; i++) {
      xv = new ScalarType[numEquations];
      bv = new ScalarType[numEquations];
      xexactv = new ScalarType[numEquations];
    }

    for (i=0; i<nrhs; i++)
      for (j=0; j<numEquations; j++) {
	xexactv[i][j] = ((ScalarType) rand())/ ((ScalarType) RAND_MAX);
	xv[i][j] = 0.0;
	bv[i][j] = 0.0;
      }
  OrdinalType * indices = 0;
  ScalarType * values = 0;
  OrdinalType * pointers = 0;
  OrdinalType * allIndices = 0;
  ScalarType * allValues = 0;
  OrdinalType * profiles = 0;

  if (generateClassicHbMatrix) {
    allIndices = new OrdinalType[numEquations*npoints];
    allValues = new ScalarType[numEquations*npoints];
    profiles = new OrdinalType[numEquations+1];
  }
  else {
  indices = new OrdinalType*[numEquations];
  values = new ScalarType*[numEquations];
  pointers = new OrdinalType[numEquations];

  OrdinalType * curIndices;
  ScalarType * curValues;

  ScalarType dnpoints = (ScalarType) npoints;
  OrdinalType numEntries = 0;
  pointers[0] = 0;
  for (i=0; i<numEquations; i++) {

    OrdinalType rowID = i;
    OrdinalType numIndices = 0;
    
    if (generateClassicHbMatrix) {
      curIndices = allIndices+numEntries;
      curValues = allValues+numEntries;
    }
    else {
      curIndices = new OrdinalType[npoints];
      curValues = new ScalarType[npoints];
      indices[i] = curIndices;
      values[i] = curValues;
    }
      

    for (j=0; j<npoints; j++) {
      OrdinalType colID = rowID + xoff[j] + nx*yoff[j]; // Compute column ID based on stencil offsets
      if (colID>-1 && colID<numEquations) {
	curIindices[numIndices] = colID;
	ScalarType value = - ((ScalarType) rand())/ ((ScalarType) RAND_MAX);
	if (colID==rowID)
	  curValues[numIndices] = dnpoints - value; // Make diagonal dominant
	else
	  curValues[numIndices] = -value;
	if (isRowOriented)
	  for (k=0; k<nrhs; k++)
	    bv[i][k] += curValues[numIndices]*xexactv[curIindices[numIndices]][k];
	else
	  for (k=0; k<nrhs; k++)
	    bv[curIindices[numIndices]][k] += curValues[numIndices]*xexactv[i][k];
	numEntries++;
	numIndices++;
      }
    }
    if (generateClassicHbMatrix) 
      pointers[i+1] = pointers[i]+numIndices;
    else
      profiles[i] = numIndices;
  }

  // Now construct Kokkos objects and register pointers
  x = new Kokkos::MultiVector<OrdinalType, ScalarType>();
  b = new Kokkos::MultiVector<OrdinalType, ScalarType>();
  xexact = new Kokkos::MultiVector<OrdinalType, ScalarType>();
  A = new Kokkos::HbMatrix<OrdinalType, ScalarType>(); // Construct matrix

  // Use the generateClassicHbMatrix as a flag to toggle how multivector
  // values are initialized also.
  if (generateClassicHbMatrix) {
    x.initializeValues(numEquations, nrhs, xv[0], numEquations, 1);
    b.initializeValues(numEquations, nrhs, bv[0], numEquations, 1);
    xexact.initializeValues(numEquations, nrhs, xexactv[0], numEquations, 1);

    A.initializeStructure(numEquations, numEquations, isRowOriented, pointers, allIndices);
    A.initializeValues(allValues);
  }
  else {
    x.initializeValues(numEquations, nrhs, xv);
    b.initializeValues(numEquations, nrhs, bv);
    xexact.initializeValues(numEquations, nrhs, xexactv);

    A.initializeStructure(numEquations, numEquations, isRowOriented, profiles, indices);
    A.initializeValues(values);
  }

  
  return;
}
