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

#include "LOCA_Bifurcation_PitchforkBordVector.H"  // Class definition

LOCA::Bifurcation::PitchforkBordVector::PitchforkBordVector(
			                const NOX::Abstract::Vector& xVec,
					const NOX::Abstract::Vector& nullVec,
					double slackVar, double bifParam) :
  LOCA::MultiVector(2,2)
{
  setVector(0, xVec);
  setVector(1, nullVec);
  setScalar(0, slackVar);
  setScalar(1, bifParam);
}

LOCA::Bifurcation::PitchforkBordVector::PitchforkBordVector(
                         const LOCA::Bifurcation::PitchforkBordVector& source, 
			 NOX::CopyType type) :
  LOCA::MultiVector(source, type)
{
}


LOCA::Bifurcation::PitchforkBordVector::~PitchforkBordVector()
{
}

NOX::Abstract::Vector& 
LOCA::Bifurcation::PitchforkBordVector::operator=(
                                                const NOX::Abstract::Vector& y)
{
  return operator=(dynamic_cast<const LOCA::Bifurcation::PitchforkBordVector&>(y));
}

LOCA::MultiVector& 
LOCA::Bifurcation::PitchforkBordVector::operator=(const LOCA::MultiVector& y)
{
  return operator=(dynamic_cast<const LOCA::Bifurcation::PitchforkBordVector&>(y));
}

LOCA::Bifurcation::PitchforkBordVector& 
LOCA::Bifurcation::PitchforkBordVector::operator=(
                              const LOCA::Bifurcation::PitchforkBordVector& y)
{ 
  LOCA::MultiVector::operator=(y);
  return *this;
}

NOX::Abstract::Vector* 
LOCA::Bifurcation::PitchforkBordVector::clone(NOX::CopyType type) const
{
  return new LOCA::Bifurcation::PitchforkBordVector(*this, type);
}

void 
LOCA::Bifurcation::PitchforkBordVector::setVec(
                                         const NOX::Abstract::Vector& xVec,
					 const NOX::Abstract::Vector& nullVec,
					 double slackVar, double bifPar)
{
  setVector(0, xVec);
  setVector(1, nullVec);
  setScalar(0, slackVar);
  setScalar(1, bifPar);
}

const NOX::Abstract::Vector& 
LOCA::Bifurcation::PitchforkBordVector::getXVec() const
{
  return getVector(0);
}

const NOX::Abstract::Vector& 
LOCA::Bifurcation::PitchforkBordVector::getNullVec() const
{
  return getVector(1);
}

double 
LOCA::Bifurcation::PitchforkBordVector::getSlackVar() const
{
  return getScalar(0);
}

double 
LOCA::Bifurcation::PitchforkBordVector::getBifParam() const
{
  return getScalar(1);
}

NOX::Abstract::Vector& 
LOCA::Bifurcation::PitchforkBordVector::getXVec()
{
  return getVector(0);
}

NOX::Abstract::Vector& 
LOCA::Bifurcation::PitchforkBordVector::getNullVec()
{
  return getVector(1);
}

double& 
LOCA::Bifurcation::PitchforkBordVector::getSlackVar()
{
  return getScalar(0);
}

double& 
LOCA::Bifurcation::PitchforkBordVector::getBifParam()
{
  return getScalar(1);
}
