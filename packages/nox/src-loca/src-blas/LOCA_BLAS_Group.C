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

#include "LOCA_BLAS_Group.H"	// class definition
#include "LOCA_Parameter_Vector.H"

using namespace LOCA;
using namespace LOCA::BLAS;

Group::Group(Interface& interface) : 
  NOX::BLAS::Group(interface), locaProblemInterface(interface), params(),
  tangentVec(interface.getInitialGuess())
{}

Group::Group(const Group& source, NOX::CopyType type) : 
  NOX::BLAS::Group(source,type), 
  locaProblemInterface(source.locaProblemInterface), params(source.params),
  tangentVec(source.tangentVec)
{}

Group::~Group() 
{}

NOX::Abstract::Group& 
LOCA::BLAS::Group::operator=(const NOX::Abstract::Group& source) {
  return operator=(dynamic_cast<const Group&>(source));
}

Abstract::Group& 
LOCA::BLAS::Group::operator=(const Abstract::Group& source) {
  return operator=(dynamic_cast<const Group&>(source));
}

Group& 
LOCA::BLAS::Group::operator=(const Group& source) {
  NOX::BLAS::Group::operator=(source);
  params = source.params;
  tangentVec = source.tangentVec;
  return *this;
}

NOX::Abstract::Group*
LOCA::BLAS::Group::clone(NOX::CopyType type) const {
  return new Group(*this, type);
}

bool
LOCA::BLAS::Group::computeF() {
  locaProblemInterface.setParams(params);
  return NOX::BLAS::Group::computeF();
}

bool
LOCA::BLAS::Group::computeJacobian() {
  locaProblemInterface.setParams(params);
  return NOX::BLAS::Group::computeJacobian();
}

bool
LOCA::BLAS::Group::computeTangent(NOX::Parameter::List& params,
                                      int paramID)
{
  bool res = computeJacobian();

  NOX::Abstract::Vector* dfdpVec = tangentVec.clone(NOX::ShapeCopy);

  res = res && computeDfDp(paramID, *dfdpVec);

  res = res && applyJacobianInverse(params, *dfdpVec, tangentVec);

  delete dfdpVec;

  return res;
}

bool 
LOCA::BLAS::Group::setParams(const ParameterVector& p) 
{
  resetIsValid();
  params = p;
  return true;
}

bool 
LOCA::BLAS::Group::computeParams(const ParameterVector& oldParams, 
				 const ParameterVector& direction, 
				 double step) 
{
  resetIsValid();
  params = oldParams;
  params.update(step, direction, 1.0);
  return true;
}

const ParameterVector& 
LOCA::BLAS::Group::getParams() const
{
  return params;
}

const NOX::Abstract::Vector&
LOCA::BLAS::Group::getTangent() const
{
  return tangentVec;
}

bool 
LOCA::BLAS::Group::print() const
{
  cout << "p = " << params << "\n";
  NOX::BLAS::Group::print();
  return true;
}
