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

#include "LOCA_Continuation_ExtendedVector.H"  // Class definition

LOCA::Continuation::ExtendedVector::ExtendedVector(
					   const NOX::Abstract::Vector& xVec,
					   double param) :
  LOCA::Extended::Vector(1,1)
{
  LOCA::Extended::Vector::setVector(0, xVec);
  LOCA::Extended::Vector::setScalar(0, param);
}

LOCA::Continuation::ExtendedVector::ExtendedVector(
			    const LOCA::Continuation::ExtendedVector& source, 
			    NOX::CopyType type) :
  LOCA::Extended::Vector(source, type)
{
}

LOCA::Continuation::ExtendedVector::~ExtendedVector()
{
}

LOCA::Extended::Vector& 
LOCA::Continuation::ExtendedVector::operator=(
					     const LOCA::Extended::Vector& y)
{
  return 
    operator=(dynamic_cast<const LOCA::Continuation::ExtendedVector&>(y));
}

NOX::Abstract::Vector& 
LOCA::Continuation::ExtendedVector::operator=(const NOX::Abstract::Vector& y)
{
 return operator=(dynamic_cast<const LOCA::Continuation::ExtendedVector&>(y));
}

LOCA::Continuation::ExtendedVector& 
LOCA::Continuation::ExtendedVector::operator=(const 
				      LOCA::Continuation::ExtendedVector& y)
{
  LOCA::Extended::Vector::operator=(y);
  return *this;
}

NOX::Abstract::Vector* 
LOCA::Continuation::ExtendedVector::clone(NOX::CopyType type) const
{
  return new LOCA::Continuation::ExtendedVector(*this, type);
}

void 
LOCA::Continuation::ExtendedVector::setVec(const NOX::Abstract::Vector& xVec,
				   double param)
{
  LOCA::Extended::Vector::setVector(0, xVec);
  LOCA::Extended::Vector::setScalar(0, param);
}

const NOX::Abstract::Vector& 
LOCA::Continuation::ExtendedVector::getXVec() const
{
  return LOCA::Extended::Vector::getVector(0);
}

double 
LOCA::Continuation::ExtendedVector::getParam() const
{
  return LOCA::Extended::Vector::getScalar(0);
}

NOX::Abstract::Vector& 
LOCA::Continuation::ExtendedVector::getXVec()
{
  return LOCA::Extended::Vector::getVector(0);
}

double& 
LOCA::Continuation::ExtendedVector::getParam()
{
  return LOCA::Extended::Vector::getScalar(0);
}
