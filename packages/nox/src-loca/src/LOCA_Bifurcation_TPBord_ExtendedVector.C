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

#include "LOCA_Bifurcation_TPBord_ExtendedVector.H"  // Class definition

LOCA::Bifurcation::TPBord::ExtendedVector::ExtendedVector(
					  const NOX::Abstract::Vector& xVec,
					  const NOX::Abstract::Vector& nullVec,
					  double bifParam) :
  LOCA::Extended::Vector(2,1)
{
  setVector(0, xVec);
  setVector(1, nullVec);
  setScalar(0, bifParam);
}

LOCA::Bifurcation::TPBord::ExtendedVector::ExtendedVector(
                      const LOCA::Bifurcation::TPBord::ExtendedVector& source,
		      NOX::CopyType type) :
  LOCA::Extended::Vector(source, type)
{
}


LOCA::Bifurcation::TPBord::ExtendedVector::~ExtendedVector()
{
}

NOX::Abstract::Vector& 
LOCA::Bifurcation::TPBord::ExtendedVector::operator=(
					      const NOX::Abstract::Vector& y)
{
  return operator=(dynamic_cast<const LOCA::Bifurcation::TPBord::ExtendedVector&>(y));
}

LOCA::Extended::Vector& 
LOCA::Bifurcation::TPBord::ExtendedVector::operator=(
					     const LOCA::Extended::Vector& y)
{
  return 
    operator=(dynamic_cast<const LOCA::Bifurcation::TPBord::ExtendedVector&>(y));
}

LOCA::Bifurcation::TPBord::ExtendedVector& 
LOCA::Bifurcation::TPBord::ExtendedVector::operator=(
                           const LOCA::Bifurcation::TPBord::ExtendedVector& y)
{ 
  LOCA::Extended::Vector::operator=(y);
  return *this;
}

NOX::Abstract::Vector* 
LOCA::Bifurcation::TPBord::ExtendedVector::clone(NOX::CopyType type) const
{
  return new LOCA::Bifurcation::TPBord::ExtendedVector(*this, type);
}

void 
LOCA::Bifurcation::TPBord::ExtendedVector::setVec(
					const NOX::Abstract::Vector& xVec,
					const NOX::Abstract::Vector& nullVec,
					double bifPar)
{
  setVector(0, xVec);
  setVector(1, nullVec);
  setScalar(0, bifPar);
}

const NOX::Abstract::Vector& 
LOCA::Bifurcation::TPBord::ExtendedVector::getXVec() const
{
  return getVector(0);
}

const NOX::Abstract::Vector& 
LOCA::Bifurcation::TPBord::ExtendedVector::getNullVec() const
{
  return getVector(1);
}

double 
LOCA::Bifurcation::TPBord::ExtendedVector::getBifParam() const
{
  return getScalar(0);
}

NOX::Abstract::Vector& 
LOCA::Bifurcation::TPBord::ExtendedVector::getXVec()
{
  return getVector(0);
}

NOX::Abstract::Vector& 
LOCA::Bifurcation::TPBord::ExtendedVector::getNullVec()
{
  return getVector(1);
}

double& 
LOCA::Bifurcation::TPBord::ExtendedVector::getBifParam()
{
  return getScalar(0);
}
