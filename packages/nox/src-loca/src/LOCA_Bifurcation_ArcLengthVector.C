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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "LOCA_Bifurcation_ArcLengthVector.H"  // Class definition

using namespace LOCA;
using namespace LOCA::Bifurcation; 

ArcLengthVector::ArcLengthVector(const NOX::Abstract::Vector& xVec,
				 double arcParam) :
  MultiVector(1,1)
{
  setVector(0, xVec);
  setScalar(0, arcParam);
}

ArcLengthVector::ArcLengthVector(const ArcLengthVector& source, 
				 NOX::CopyType type) :
  MultiVector(source, type)
{
}

ArcLengthVector::~ArcLengthVector()
{
}

NOX::Abstract::Vector& ArcLengthVector::operator=(const NOX::Abstract::Vector& y)
{
 return operator=(dynamic_cast<const ArcLengthVector&>(y));
}

ArcLengthVector& ArcLengthVector::operator=(const ArcLengthVector& y)
{
  MultiVector::operator=(y);
  return *this;
}

NOX::Abstract::Vector* ArcLengthVector::clone(NOX::CopyType type) const
{
  return new ArcLengthVector(*this, type);
}

void ArcLengthVector::setVec(const NOX::Abstract::Vector& xVec,
                             double arcPar)
{
  setVector(0, xVec);
  setScalar(0, arcPar);
}

const NOX::Abstract::Vector& ArcLengthVector::getXVec() const
{
  return getVector(0);
}

double ArcLengthVector::getArcParam() const
{
  return getScalar(0);
}

NOX::Abstract::Vector& ArcLengthVector::getXVec()
{
  return getVector(0);
}

double& ArcLengthVector::getArcParam()
{
  return getScalar(0);
}
