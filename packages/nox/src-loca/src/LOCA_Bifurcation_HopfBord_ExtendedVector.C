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

#include "LOCA_Bifurcation_HopfBord_ExtendedVector.H"  // Class definition

LOCA::Bifurcation::HopfBord::ExtendedVector::ExtendedVector(
			       const NOX::Abstract::Vector& xVec,
			       const NOX::Abstract::Vector& realEigenVec,
		               const NOX::Abstract::Vector& imagEigenVec,
			       double frequency,
			       double bifParam) :
  LOCA::Extended::Vector(3,2)
{
  setVector(0, xVec);
  setVector(1, realEigenVec);
  setVector(2, imagEigenVec);
  setScalar(0, frequency);
  setScalar(1, bifParam);
}

LOCA::Bifurcation::HopfBord::ExtendedVector::ExtendedVector(
                    const LOCA::Bifurcation::HopfBord::ExtendedVector& source,
		    NOX::CopyType type) :
  LOCA::Extended::Vector(source, type)
{
}


LOCA::Bifurcation::HopfBord::ExtendedVector::~ExtendedVector()
{
}

NOX::Abstract::Vector& 
LOCA::Bifurcation::HopfBord::ExtendedVector::operator=(
					      const NOX::Abstract::Vector& y)
{
  return operator=(dynamic_cast<const LOCA::Bifurcation::HopfBord::ExtendedVector&>(y));
}

LOCA::Extended::Vector& 
LOCA::Bifurcation::HopfBord::ExtendedVector::operator=(
					      const LOCA::Extended::Vector& y)
{
  return operator=(dynamic_cast<const LOCA::Bifurcation::HopfBord::ExtendedVector&>(y));
}

LOCA::Bifurcation::HopfBord::ExtendedVector& 
LOCA::Bifurcation::HopfBord::ExtendedVector::operator=(
                         const LOCA::Bifurcation::HopfBord::ExtendedVector& y)
{ 
  LOCA::Extended::Vector::operator=(y);
  return *this;
}

NOX::Abstract::Vector* 
LOCA::Bifurcation::HopfBord::ExtendedVector::clone(NOX::CopyType type) const
{
  return new LOCA::Bifurcation::HopfBord::ExtendedVector(*this, type);
}

void 
LOCA::Bifurcation::HopfBord::ExtendedVector::setVec(
			       const NOX::Abstract::Vector& xVec,
			       const NOX::Abstract::Vector& realEigenVec,
			       const NOX::Abstract::Vector& imagEigenVec,
			       double frequency,
			       double bifPar)
{
  setVector(0, xVec);
  setVector(1, realEigenVec);
  setVector(2, imagEigenVec);
  setScalar(0, frequency);
  setScalar(1, bifPar);
}

const NOX::Abstract::Vector& 
LOCA::Bifurcation::HopfBord::ExtendedVector::getXVec() const
{
  return getVector(0);
}

const NOX::Abstract::Vector& 
LOCA::Bifurcation::HopfBord::ExtendedVector::getRealEigenVec() const
{
  return getVector(1);
}

const NOX::Abstract::Vector& 
LOCA::Bifurcation::HopfBord::ExtendedVector::getImagEigenVec() const
{
  return getVector(2);
}

double 
LOCA::Bifurcation::HopfBord::ExtendedVector::getFrequency() const
{
  return getScalar(0);
}

double 
LOCA::Bifurcation::HopfBord::ExtendedVector::getBifParam() const
{
  return getScalar(1);
}

NOX::Abstract::Vector& 
LOCA::Bifurcation::HopfBord::ExtendedVector::getXVec()
{
  return getVector(0);
}

NOX::Abstract::Vector& 
LOCA::Bifurcation::HopfBord::ExtendedVector::getRealEigenVec()
{
  return getVector(1);
}

NOX::Abstract::Vector& 
LOCA::Bifurcation::HopfBord::ExtendedVector::getImagEigenVec()
{
  return getVector(2);
}

double& 
LOCA::Bifurcation::HopfBord::ExtendedVector::getFrequency()
{
  return getScalar(0);
}

double& 
LOCA::Bifurcation::HopfBord::ExtendedVector::getBifParam()
{
  return getScalar(1);
}
