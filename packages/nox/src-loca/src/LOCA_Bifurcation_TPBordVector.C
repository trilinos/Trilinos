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

#include "LOCA_Bifurcation_TPBordVector.H"  // Class definition

using namespace LOCA;
using namespace LOCA::Bifurcation; 

TPBordVector::TPBordVector(const NOX::Abstract::Vector& xVec,
			   const NOX::Abstract::Vector& nullVec,
			   double bifParam) :
  MultiVector(2,1)
{
  setVector(0, xVec);
  setVector(1, nullVec);
  setScalar(0, bifParam);
}

TPBordVector::TPBordVector(const TPBordVector& source, NOX::CopyType type) :
  MultiVector(source, type)
{
}


TPBordVector::~TPBordVector()
{
}

NOX::Abstract::Vector& TPBordVector::operator=(const NOX::Abstract::Vector& y)
{
  return operator=(dynamic_cast<const TPBordVector&>(y));
}

LOCA::MultiVector& TPBordVector::operator=(const LOCA::MultiVector& y)
{
  return operator=(dynamic_cast<const TPBordVector&>(y));
}

TPBordVector& TPBordVector::operator=(const TPBordVector& y)
{ 
  MultiVector::operator=(y);
  return *this;
}

NOX::Abstract::Vector* TPBordVector::clone(NOX::CopyType type) const
{
  return new TPBordVector(*this, type);
}

void TPBordVector::setVec(const NOX::Abstract::Vector& xVec,
			  const NOX::Abstract::Vector& nullVec,
			  double bifPar)
{
  setVector(0, xVec);
  setVector(1, nullVec);
  setScalar(0, bifPar);
}

const NOX::Abstract::Vector& TPBordVector::getXVec() const
{
  return getVector(0);
}

const NOX::Abstract::Vector& TPBordVector::getNullVec() const
{
  return getVector(1);
}

double TPBordVector::getBifParam() const
{
  return getScalar(0);
}

NOX::Abstract::Vector& TPBordVector::getXVec()
{
  return getVector(0);
}

NOX::Abstract::Vector& TPBordVector::getNullVec()
{
  return getVector(1);
}

double& TPBordVector::getBifParam()
{
  return getScalar(0);
}
