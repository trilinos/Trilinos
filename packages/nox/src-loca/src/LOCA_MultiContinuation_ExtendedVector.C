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

#include "LOCA_MultiContinuation_ExtendedVector.H"  // Class definition
#include "LOCA_MultiContinuation_ExtendedMultiVector.H"

LOCA::MultiContinuation::ExtendedVector::ExtendedVector(
					   const NOX::Abstract::Vector& xVec,
					   int nScalars) :
  LOCA::Extended::Vector(1,nScalars)
{
  LOCA::Extended::Vector::setVector(0, xVec);
}

LOCA::MultiContinuation::ExtendedVector::ExtendedVector(
			const LOCA::MultiContinuation::ExtendedVector& source, 
			NOX::CopyType type) :
  LOCA::Extended::Vector(source, type)
{
}

LOCA::MultiContinuation::ExtendedVector::~ExtendedVector()
{
}

LOCA::Extended::Vector& 
LOCA::MultiContinuation::ExtendedVector::operator=(
					     const LOCA::Extended::Vector& y)
{
  return 
    operator=(dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(y));
}

NOX::Abstract::Vector& 
LOCA::MultiContinuation::ExtendedVector::operator=(
					       const NOX::Abstract::Vector& y)
{
 return operator=(dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(y));
}

LOCA::MultiContinuation::ExtendedVector& 
LOCA::MultiContinuation::ExtendedVector::operator=(const 
				   LOCA::MultiContinuation::ExtendedVector& y)
{
  LOCA::Extended::Vector::operator=(y);
  return *this;
}

NOX::Abstract::Vector* 
LOCA::MultiContinuation::ExtendedVector::clone(NOX::CopyType type) const
{
  return new LOCA::MultiContinuation::ExtendedVector(*this, type);
}

const NOX::Abstract::Vector& 
LOCA::MultiContinuation::ExtendedVector::getXVec() const
{
  return LOCA::Extended::Vector::getVector(0);
}

NOX::Abstract::Vector& 
LOCA::MultiContinuation::ExtendedVector::getXVec()
{
  return LOCA::Extended::Vector::getVector(0);
}

LOCA::MultiContinuation::ExtendedVector::ExtendedVector(int nScalars) :
  LOCA::Extended::Vector(1,nScalars)
{
}

LOCA::Extended::MultiVector*
LOCA::MultiContinuation::ExtendedVector::generateMultiVector(int nColumns, 
							int nVectorRows, 
							int nScalarRows) const
{
  return new LOCA::MultiContinuation::ExtendedMultiVector(nColumns,
							  nScalarRows);
}
