// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "NOX_Abstract_Vector.H"

// Included multivector declarations if requested
#ifdef HAVE_NOX_MULTIVECS
#include "NOX_MultiVector.H"
#endif

NOX::Abstract::Vector& NOX::Abstract::Vector::random(bool useSeed, int seed) 
{
  cerr << "NOX::Abstract::Vector::random() function not implemented" << endl;
  throw "NOX Error";
  return *this;
}

void NOX::Abstract::Vector::print() const
{
  return;
}

#ifdef HAVE_NOX_MULTIVECS

NOX::Abstract::MultiVector* 
NOX::Abstract::Vector::createMultiVector(
				    const NOX::Abstract::Vector* const* vecs,
				    int numVecs, NOX::CopyType type) const
{
  if (numVecs < 0) {
    cerr << "NOX::Abstract::Vector::createMultiVector:  Error!  Multivector" 
	 << " must have postive number of columns!" << endl;
    throw "NOX Error";
  }

  const NOX::Abstract::Vector** tmp = 
    new const NOX::Abstract::Vector*[numVecs+1];

  tmp[0] = this;
  for (int i=0; i<numVecs; i++)
    tmp[i+1] = vecs[i];

  NOX::MultiVector* mv = new NOX::MultiVector(tmp, numVecs+1, type);

  delete [] tmp;

  return mv;
}

#endif
