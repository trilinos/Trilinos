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

#include "LOCA_LAPACK_Group.H"	// class definition

LOCA::LAPACK::Group::Group(Interface& interface) : 
  NOX::LAPACK::Group(interface), 
  locaProblemInterface(interface), 
  params(),
  tangentVec(interface.getInitialGuess()),
  isValidTangent(false)
{}

LOCA::LAPACK::Group::Group::Group(const Group& source, NOX::CopyType type) : 
  NOX::LAPACK::Group(source,type), 
  locaProblemInterface(source.locaProblemInterface), 
  params(source.params),
  tangentVec(source.tangentVec,type)
{
  switch (type) {
    
  case NOX::DeepCopy:
    
    isValidTangent = source.isValidTangent;
    break;

  case NOX::ShapeCopy:
    resetIsValid();
    break;

  default:
    cerr << "LOCA::LAPACK::Group - invalid CopyType for copy constructor." << endl;
    throw "LOCA LAPACK Error";
  }

}

LOCA::LAPACK::Group::~Group() 
{}

NOX::Abstract::Group& 
LOCA::LAPACK::Group::operator=(const NOX::Abstract::Group& source) {
  return operator=(dynamic_cast<const LOCA::LAPACK::Group&>(source));
}

LOCA::Abstract::Group& 
LOCA::LAPACK::Group::operator=(const LOCA::Abstract::Group& source) {
  return operator=(dynamic_cast<const LOCA::LAPACK::Group&>(source));
}

NOX::LAPACK::Group&
LOCA::LAPACK::Group::operator=(const NOX::LAPACK::Group& source) {
  return operator=(dynamic_cast<const LOCA::LAPACK::Group&>(source));
}

LOCA::LAPACK::Group& 
LOCA::LAPACK::Group::operator=(const LOCA::LAPACK::Group& source) {
  NOX::LAPACK::Group::operator=(source);
  params = source.params;
  tangentVec = source.tangentVec;
  isValidTangent = source.isValidTangent;
  return *this;
}

NOX::Abstract::Group*
LOCA::LAPACK::Group::clone(NOX::CopyType type) const {
  return new LOCA::LAPACK::Group(*this, type);
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::computeF() {
  locaProblemInterface.setParams(params);
  return NOX::LAPACK::Group::computeF();
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::computeJacobian() {
  locaProblemInterface.setParams(params);
  return NOX::LAPACK::Group::computeJacobian();
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::computeTangent(NOX::Parameter::List& params,
                                      int paramID)
{
  if (isValidTangent)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType res = computeJacobian();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  NOX::Abstract::Vector* dfdpVec = tangentVec.clone(NOX::ShapeCopy);

  res = computeDfDp(paramID, *dfdpVec);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  res = applyJacobianInverse(params, *dfdpVec, tangentVec);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  tangentVec.scale(-1.0);

  delete dfdpVec;

  return res;
}

void
LOCA::LAPACK::Group::setParams(const LOCA::ParameterVector& p) 
{
  resetIsValid();
  params = p;
}

void
LOCA::LAPACK::Group::computeParams(const LOCA::ParameterVector& oldParams, 
				   const LOCA::ParameterVector& direction, 
				   double step) 
{
  resetIsValid();
  params = oldParams;
  params.update(step, direction, 1.0);
}

const LOCA::ParameterVector& 
LOCA::LAPACK::Group::getParams() const
{
  return params;
}

const NOX::Abstract::Vector&
LOCA::LAPACK::Group::getTangent() const
{
  return tangentVec;
}

void 
LOCA::LAPACK::Group::print() const
{
  cout << "p = " << params << "\n";
  NOX::LAPACK::Group::print();
}

void LOCA::LAPACK::Group::resetIsValid()
{
  isValidTangent = false;
  NOX::LAPACK::Group::resetIsValid();
}
