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

#include "LOCA_Abstract_Group.H"
#include "LOCA_Parameter_Vector.H"

using namespace LOCA;
using namespace Abstract;
using NOX::Abstract::Vector;

LOCA::Abstract::Group::Group(const DerivUtils& d)
  : derivPtr(d.clone(NOX::DeepCopy))
{
}

LOCA::Abstract::Group::Group(const Group& source, NOX::CopyType type)
  : derivPtr(source.derivPtr->clone(type))
{
}


LOCA::Abstract::Group::~Group() 
{
  delete derivPtr;
}

LOCA::Abstract::Group&
LOCA::Abstract::Group::operator=(const LOCA::Abstract::Group& source)
{

  // Protect against A = A
  if (this != &source) {
    NOX::CopyType type = NOX::DeepCopy;

    // Delete old values
    delete derivPtr;

    // Copy values
    derivPtr = source.derivPtr->clone(type);
  }

  return *this;
}

bool
LOCA::Abstract::Group::computeTangent(NOX::Parameter::List& params, 
                                      int paramID, 
                                      NOX::Abstract::Vector& result) 
{
  bool res = computeJacobian();

  NOX::Abstract::Vector *dfdpVec = result.clone(NOX::ShapeCopy);

  res = res && computeDfDp(paramID, *dfdpVec);

  res = res && applyJacobianInverse(params, *dfdpVec, result);

  delete dfdpVec;

  return res;
}

bool
LOCA::Abstract::Group::computeDfDp(int paramID, NOX::Abstract::Vector& result) 
{
  return derivPtr->computeDfDp(*this,paramID, result);
}

//NOX::Abstract::Group*
//LOCA::Abstract::Group::clone(NOX::CopyType type) const 
//{
//  return new Group(*this, type);
//}
