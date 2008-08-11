//@HEADER
// ************************************************************************
// 
//                 Kokkos: A Fast Kernel Package
//              Copyright (2004) Sandia Corporation
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

#ifndef KOKKOS_TEST_GENERATEOSKIPROBLEM_H
#define KOKKOS_TEST_GENERATEOSKIPROBLEM_H
#include "Kokkos_OskiVector.hpp"
#include "Kokkos_OskiMultiVector.hpp"
#include "Kokkos_OskiMatrix.hpp"

namespace KokkosTest {


//! KokkosTest::GenerateOskiProblem: Constructs a 2D PDE finite difference matrix using the list of x and y offsets.
  /*!
 
  nx      (In) - number of grid points in x direction

  ny      (In) - number of grid points in y direction

  The total number of equations will be nx*ny ordered such that the x direction changes
  most rapidly: 
  First equation is at point (0,0)
  Second at                  (1,0)
  ...
  nx equation at             (nx-1,0)
  nx+1st equation at         (0,1)

  npoints (In) - number of points in finite difference stencil
  xoff    (In) - stencil offsets in x direction (of length npoints)
  yoff    (In) - stencil offsets in y direction (of length npoints)
  A standard 5-point finite difference stencil would be described as:
  npoints = 5
  xoff = [-1, 1, 0,  0, 0]
  yoff = [ 0, 0, 0, -1, 1]

  nrhs - Number of rhs to generate. (First interface produces vectors, so nrhs is not needed

  A      (Out) - Kokkos::OskiMatrix constructed for nx by ny grid using prescribed stencil
  Off-diagonal values are random between 0 and 1.  If diagonal is part of stencil,
  diagonal will be slightly diag dominant.
  x      (Out) - Initial guess vector set to zero.
  b      (Out) - Generated RHS.  Values satisfy b = A*xexact
  xexact (Out) - Generated exact solution to Ax = b.

  Note: Caller of this function is responsible for deleting all output objects.
  */

  template<typename OrdinalType, typename ScalarType>
  class GenerateOskiProblem {
  public:

    //@{ \name Constructors/Destructor.

    //! Single RHS constuctor.
/*    GenerateOskiProblem(bool generateClassicHbMatrix, bool isRowOriented, bool hasImplicitUnitDiagonal,
		      OrdinalType nx, OrdinalType ny, OrdinalType npoints, 
		      OrdinalType * xoff, OrdinalType * yoff,
		      Kokkos::CisMatrix<OrdinalType, ScalarType> *& A, 
		      Kokkos::Vector<OrdinalType, ScalarType> *& x, 
		      Kokkos::Vector<OrdinalType, ScalarType> *& b,
		      Kokkos::Vector<OrdinalType, ScalarType> *&xexact,
		      OrdinalType & numEntries);
*/
    //! Multi RHS constuctor.
    GenerateOskiProblem(bool isRowOriented, bool hasImplicitUnitDiagonal,
		      OrdinalType nx, OrdinalType ny, OrdinalType npoints, 
		      OrdinalType * xoff, OrdinalType * yoff, OrdinalType nrhs,
		      Kokkos::CisMatrix<OrdinalType, ScalarType> *& A, 
		      Kokkos::MultiVector<OrdinalType, ScalarType> *& x, 
		      Kokkos::MultiVector<OrdinalType, ScalarType> *& b,
		      Kokkos::MultiVector<OrdinalType, ScalarType> *&xexact,
		      OrdinalType & numEntries);

  public:

    //! HbMatrix Destructor
    virtual ~GenerateOskiProblem();
    //@}
  
   private:

    //! Copy constructor (Not implemented).
    GenerateOskiProblem(const GenerateOskiProblem<OrdinalType, ScalarType>& source){};
/*
   void GenerateProblem(bool generateClassicHbMatrix, bool isRowOriented, bool hasImplicitUnitDiagonal,
			   OrdinalType nx, OrdinalType ny, OrdinalType npoints, 
			   OrdinalType * xoff, OrdinalType * yoff,
			   Kokkos::CisMatrix<OrdinalType, ScalarType> *& A, 
			   Kokkos::Vector<OrdinalType, ScalarType> *& x, 
			   Kokkos::Vector<OrdinalType, ScalarType> *& b,
			   Kokkos::Vector<OrdinalType, ScalarType> *&xexact,
			   OrdinalType & numEntries);
*/
      void GenerateProblem(bool isRowOriented, bool hasImplicitUnitDiagonal,
		       OrdinalType nx, OrdinalType ny, OrdinalType npoints, 
		       OrdinalType * xoff, OrdinalType * yoff, OrdinalType nrhs,
		       Kokkos::CisMatrix<OrdinalType, ScalarType> *& A, 
		       Kokkos::MultiVector<OrdinalType, ScalarType> *& x, 
		       Kokkos::MultiVector<OrdinalType, ScalarType> *& b,
		       Kokkos::MultiVector<OrdinalType, ScalarType> *&xexact,
			     OrdinalType & numEntries);

  private:

    Kokkos::OskiVector<OrdinalType, ScalarType> * xd1;
    Kokkos::OskiVector<OrdinalType, ScalarType> *bd1;
    Kokkos::OskiVector<OrdinalType, ScalarType> *xexactd1;
    Kokkos::OskiMultiVector<OrdinalType, ScalarType> * xd;
    Kokkos::OskiMultiVector<OrdinalType, ScalarType> * bd;
    Kokkos::OskiMultiVector<OrdinalType, ScalarType> * xexactd;
    Kokkos::OskiMatrix<OrdinalType, ScalarType> * Ad; // Construct matrix
    ScalarType ** xv;
    ScalarType ** bv;
    ScalarType ** xexactv;
    OrdinalType numEquations;
    OrdinalType nrhsv;

    OrdinalType ** indices;
    ScalarType ** values;
    OrdinalType * pointers;
    OrdinalType * allIndices;
    ScalarType * allValues;
    OrdinalType * profiles;
    //bool generateClassicHbMatrix_;
    bool hasImplicitUnitDiagonal_;
  };
} // namespace KokkosTest

using namespace KokkosTest;
//==============================================================================
/*template<typename OrdinalType, typename ScalarType>
GenerateOskiProblem<OrdinalType, ScalarType>::
GenerateOskiProblem(bool generateClassicHbMatrix, bool isRowOriented, bool hasImplicitUnitDiagonal,
		      OrdinalType nx, OrdinalType ny, OrdinalType npoints, 
		      OrdinalType * xoff, OrdinalType * yoff,
		      Kokkos::CisMatrix<OrdinalType, ScalarType> *& A, 
		      Kokkos::Vector<OrdinalType, ScalarType> *& x, 
		      Kokkos::Vector<OrdinalType, ScalarType> *& b,
		      Kokkos::Vector<OrdinalType, ScalarType> *&xexact,
		      OrdinalType & numEntries):
  xd1(0),
  bd1(0),
  xexactd1(0),
  xd(0),
  bd(0),
  xexactd(0),
  Ad(0),
  xv(0),
  bv(0),
  xexactv(0),
  numEquations(0),
  nrhsv(0),
  indices(0),
  values(0),
  pointers(0),
  allIndices(0),
  allValues(0),
  profiles(0),
  generateClassicHbMatrix_(generateClassicHbMatrix),
  hasImplicitUnitDiagonal_(hasImplicitUnitDiagonal)
 {
  GenerateProblem(generateClassicHbMatrix, isRowOriented, hasImplicitUnitDiagonal,
		  nx, ny, npoints, xoff, yoff, A, x, b, xexact, numEntries);
}
  */
//==============================================================================
template<typename OrdinalType, typename ScalarType>
GenerateOskiProblem<OrdinalType, ScalarType>::
GenerateOskiProblem(bool isRowOriented, bool hasImplicitUnitDiagonal,
		      OrdinalType nx, OrdinalType ny, OrdinalType npoints, 
		      OrdinalType * xoff, OrdinalType * yoff, OrdinalType nrhs,
		      Kokkos::CisMatrix<OrdinalType, ScalarType> *& A, 
		      Kokkos::MultiVector<OrdinalType, ScalarType> *& x, 
		      Kokkos::MultiVector<OrdinalType, ScalarType> *& b,
		      Kokkos::MultiVector<OrdinalType, ScalarType> *&xexact,
		      OrdinalType & numEntries):
  xd1(0),
  bd1(0),
  xexactd1(0),
  xd(0),
  bd(0),
  xexactd(0),
  Ad(0),
  xv(0),
  bv(0),
  xexactv(0),
  numEquations(0),
  nrhsv(0),
  indices(0),
  values(0),
  pointers(0),
  allIndices(0),
  allValues(0),
  profiles(0),
  //generateClassicHbMatrix_(generateClassicHbMatrix),
  hasImplicitUnitDiagonal_(hasImplicitUnitDiagonal) {
  GenerateProblem(isRowOriented,  hasImplicitUnitDiagonal, 
		  nx, ny, npoints, xoff, yoff, nrhs, A, x, b, xexact, numEntries);
}
  
//==============================================================================
/*template<typename OrdinalType, typename ScalarType>
void GenerateOskiProblem<OrdinalType, ScalarType>::
GenerateProblem(bool generateClassicHbMatrix, bool isRowOriented, bool hasImplicitUnitDiagonal,
		OrdinalType nx, OrdinalType ny, OrdinalType npoints, 
		OrdinalType * xoff, OrdinalType * yoff,
		Kokkos::CisMatrix<OrdinalType, ScalarType> *& A, 
		Kokkos::Vector<OrdinalType, ScalarType> *& x, 
		Kokkos::Vector<OrdinalType, ScalarType> *& b,
		Kokkos::Vector<OrdinalType, ScalarType> *&xexact,
		OrdinalType & numEntries) {

  Kokkos::MultiVector<OrdinalType, ScalarType> * x1, * b1, * xexact1;
	
  GenerateProblem(generateClassicHbMatrix, isRowOriented,  hasImplicitUnitDiagonal,
		  nx, ny, npoints, xoff, yoff, 1, A, x1, b1, xexact1, numEntries);

  xd1 = new Kokkos::DenseVector<OrdinalType, ScalarType>();
  bd1 = new Kokkos::DenseVector<OrdinalType, ScalarType>();
  xexactd1 = new Kokkos::DenseVector<OrdinalType, ScalarType>();

  xd1->initializeValues(x1->getNumRows(), x1->getValues(0), x1->getColInc());
  bd1->initializeValues(b1->getNumRows(), b1->getValues(0), b1->getColInc());
  xexactd1->initializeValues(xexact1->getNumRows(), xexact1->getValues(0), xexact1->getColInc());

  x = dynamic_cast<Kokkos::Vector<OrdinalType, ScalarType> *>(xd1);
  b = dynamic_cast<Kokkos::Vector<OrdinalType, ScalarType> *>(bd1);
  xexact = dynamic_cast<Kokkos::Vector<OrdinalType, ScalarType> *>(xexactd1);

  return;
}
*/
template<typename OrdinalType, typename ScalarType>
void GenerateOskiProblem<OrdinalType, ScalarType>::
GenerateProblem(bool isRowOriented, bool hasImplicitUnitDiagonal,
		       OrdinalType nx, OrdinalType ny, OrdinalType npoints, 
		       OrdinalType * xoff, OrdinalType * yoff, OrdinalType nrhs,
		       Kokkos::CisMatrix<OrdinalType, ScalarType> *& A, 
		       Kokkos::MultiVector<OrdinalType, ScalarType> *& x, 
		       Kokkos::MultiVector<OrdinalType, ScalarType> *& b,
		       Kokkos::MultiVector<OrdinalType, ScalarType> *&xexact,
		       OrdinalType & numEntries) {

  // Number of equations is nx*ny.
  numEquations = nx*ny;
  nrhsv = nrhs;
  OrdinalType i, j, k;

  

  xv = new ScalarType*[nrhs];
  bv = new ScalarType*[nrhs];
  xexactv = new ScalarType*[nrhs];

    xv[0] = new ScalarType[numEquations*nrhs];
    bv[0] = new ScalarType[numEquations*nrhs];
    xexactv[0] = new ScalarType[numEquations*nrhs];
    for (i=1; i<nrhs; i++) {
    xv[i] = xv[0]+numEquations*i;
    bv[i] = bv[0]+numEquations*i;
    xexactv[i] = xexactv[0]+numEquations*i;
    }

  for (i=0; i<nrhs; i++)
    for (j=0; j<numEquations; j++) {
      xexactv[i][j] = ((ScalarType) rand())/ ((ScalarType) RAND_MAX);
      xv[i][j] = 0.0;
      bv[i][j] = 0.0;
    }
  indices = 0;
  values = 0;
  pointers = 0;
  allIndices = 0;
  allValues = 0;
  profiles = 0;

    allIndices = new OrdinalType[numEquations*npoints];
    allValues = new ScalarType[numEquations*npoints];
    pointers = new OrdinalType[numEquations+1];

  OrdinalType * curIndices;
  ScalarType * curValues;

  ScalarType dnpoints = (ScalarType) npoints;
  numEntries = 0;
  pointers[0] = 0;
  for (i=0; i<numEquations; i++) {

    OrdinalType rowID = i;
    OrdinalType numIndices = 0;
    
      curIndices = allIndices+numEntries;
      curValues = allValues+numEntries;
      
    if (hasImplicitUnitDiagonal_) 
      for (k=0; k<nrhs; k++)
	bv[k][i] += xexactv[k][i];
    for (j=0; j<npoints; j++) {
      OrdinalType colID = rowID + xoff[j] + nx*yoff[j]; // Compute column ID based on stencil offsets
      if (colID>-1 && colID<numEquations) {
	curIndices[numIndices] = colID;
	ScalarType value = - ((ScalarType) rand())/ ((ScalarType) RAND_MAX);
	if (hasImplicitUnitDiagonal_) value /= dnpoints; // scale off diagonal if diagonal is unit
	if (colID==rowID)
	  curValues[numIndices] = dnpoints - value; // Make diagonal dominant
	else
	  curValues[numIndices] = -value;
	if (isRowOriented)
	  for (k=0; k<nrhs; k++)
	    bv[k][i] += curValues[numIndices]*xexactv[k][curIndices[numIndices]];
	else
	  for (k=0; k<nrhs; k++)
	    bv[k][curIndices[numIndices]] += curValues[numIndices]*xexactv[k][i];
	numEntries++;
	numIndices++;
      }
    }
      pointers[i+1] = pointers[i]+numIndices;
  }

  // Now construct Kokkos objects and register pointers
  xd = new Kokkos::OskiMultiVector<OrdinalType, ScalarType>();
  bd = new Kokkos::OskiMultiVector<OrdinalType, ScalarType>();
  xexactd = new Kokkos::OskiMultiVector<OrdinalType, ScalarType>();
  Ad = new Kokkos::OskiMatrix<OrdinalType, ScalarType>(); // Construct matrix

    xd->initializeValues(numEquations, nrhs, xv[0], numEquations, 1);
    bd->initializeValues(numEquations, nrhs, bv[0], numEquations, 1);
    xexactd->initializeValues(numEquations, nrhs, xexactv[0], numEquations, 1);

    Ad->initializeStructure(numEquations, numEquations, isRowOriented, pointers, allIndices);
    Ad->initializeValues(allValues);

  x = dynamic_cast<Kokkos::MultiVector<OrdinalType, ScalarType> *>(xd);
  b = dynamic_cast<Kokkos::MultiVector<OrdinalType, ScalarType> *>(bd);
  xexact = dynamic_cast<Kokkos::MultiVector<OrdinalType, ScalarType> *>(xexactd);
  A = dynamic_cast<Kokkos::CisMatrix<OrdinalType, ScalarType> *>(Ad);
  return;
}


template<typename OrdinalType, typename ScalarType>
GenerateOskiProblem<OrdinalType, ScalarType>::
~GenerateOskiProblem() {

  OrdinalType i;
  if (xd1!=0) delete xd1;
  if (bd1!=0) delete bd1;
  if (xexactd1!=0) delete xexactd1;
  if (xd!=0) delete xd;
  if (bd!=0) delete bd;
  if (xexactd!=0) delete xexactd;
  if (Ad!=0) delete Ad;

  if (xv!=0) {
    delete [] xv[0];
    delete [] xv;
  }
  if (bv!=0) {
    delete [] bv[0];
    delete [] bv;
  }
  if (xexactv!=0) {
    delete [] xexactv[0];
    delete [] xexactv;
  }
  if (indices!=0) {
    for (i=0; i<numEquations; i++) if (indices[i]!=0) delete [] indices[i];
    delete [] indices;
  }
  if (values!=0) {
    for (i=0; i<numEquations; i++) if (values[i]!=0) delete [] values[i];
    delete [] values;
  }
  if (pointers!=0) delete [] pointers;
  if (allIndices!=0) delete [] allIndices;
  if (allValues!=0) delete [] allValues;
  if (profiles!=0) delete [] profiles;
}
#endif /* KOKKOS_TEST_GENERATEOSKIPROBLEM_H */
