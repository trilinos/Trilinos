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
#include "NOX_LAPACK_Wrappers.H"

LOCA::LAPACK::Group::Group(LOCA::LAPACK::Interface& interface) : 
  NOX::LAPACK::Group(interface), 
  locaProblemInterface(interface), 
  params(),
  scaleVec(dynamic_cast<const NOX::LAPACK::Vector&>(getX()))
{
  computeScaleVec();  // use default method for computing scale vector
}

LOCA::LAPACK::Group::Group(LOCA::LAPACK::Interface& interface, 
			   const NOX::LAPACK::Vector& s) : 
  NOX::LAPACK::Group(interface), 
  locaProblemInterface(interface), 
  params(),
  scaleVec(s)
{
}

LOCA::LAPACK::Group::Group(const LOCA::LAPACK::Group& source, 
			   NOX::CopyType type) : 
  NOX::LAPACK::Group(source,type), 
  locaProblemInterface(source.locaProblemInterface), 
  params(source.params),
  scaleVec(source.scaleVec)
{
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
  scaleVec = source.scaleVec;
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
LOCA::LAPACK::Group::applyJacobianInverseMulti(NOX::Parameter::List& params,
			    const NOX::Abstract::Vector* const* inputs,
			    NOX::Abstract::Vector** outputs, int nVecs) const
{
  if (nVecs < 1)
    return NOX::Abstract::Group::Failed;

  int n = inputs[0]->length();
  int m = nVecs;
  int info;

  // Copy all input vectors into one matrix
  NOX::LAPACK::Matrix B(n,m);
  const NOX::LAPACK::Vector* constVecPtr;
  for (int j=0; j<m; j++) {
    constVecPtr = dynamic_cast<const NOX::LAPACK::Vector*>(inputs[j]);
    for (int i=0; i<n; i++)
      B(i,j) = (*constVecPtr)(i);
  }

  // Compute Jacobian LU factorization if invalid
  if (!NOX::LAPACK::Group::isValidJacobianLUFact) {
    NOX::LAPACK::Group::jacobianLUFact = NOX::LAPACK::Group::jacobianMatrix;
    DGETRF_F77(&n, &n, &jacobianLUFact(0,0), &n, 
	       &NOX::LAPACK::Group::pivots[0], &info);

    if (info != 0)
      return NOX::Abstract::Group::Failed;

    NOX::LAPACK::Group::isValidJacobianLUFact = true;
  }

  // Backsolve using LU factorization
  DGETRS_F77("N", &n, &m, &jacobianLUFact(0,0), &n, &pivots[0], 
  	     &B(0,0), &n, &info);

  if (info != 0)
      return NOX::Abstract::Group::Failed;

  // Copy result from matrix
  NOX::LAPACK::Vector* vecPtr;
  for (int j=0; j<m; j++) {
    vecPtr = dynamic_cast<NOX::LAPACK::Vector*>(outputs[j]);
    for (int i=0; i<n; i++)
      (*vecPtr)(i) = B(i,j);
  }

  return NOX::Abstract::Group::Ok;
}

void
LOCA::LAPACK::Group::setParams(const LOCA::ParameterVector& p) 
{
  resetIsValid();
  params = p;
}

const LOCA::ParameterVector& 
LOCA::LAPACK::Group::getParams() const
{
  return params;
}

void
LOCA::LAPACK::Group::setParam(int paramID, double val)
{
  resetIsValid();
  params.setValue(paramID, val);
}

double
LOCA::LAPACK::Group::getParam(int paramID) const
{
  return params.getValue(paramID);
}

void
LOCA::LAPACK::Group::setParam(string paramID, double val)
{
  resetIsValid();
  params.setValue(paramID, val);
}

double
LOCA::LAPACK::Group::getParam(string paramID) const
{
  return params.getValue(paramID);
}

void 
LOCA::LAPACK::Group::printSolution(const double conParam) const
{
   printSolution(xVector, conParam);
}

void
LOCA::LAPACK::Group::printSolution(const NOX::LAPACK::Vector& x_,
                                   const double conParam) const
{
   locaProblemInterface.printSolution(x_, conParam);
}

void
LOCA::LAPACK::Group::printSolution(const NOX::Abstract::Vector& x_,
                                   const double conParam) const
{
   printSolution(dynamic_cast<const NOX::LAPACK::Vector&>(x_), conParam);
}

void
LOCA::LAPACK::Group::setScaleVec(const NOX::Abstract::Vector& s) {
  setScaleVec( dynamic_cast<const NOX::LAPACK::Vector&>(s) );
}

void
LOCA::LAPACK::Group::setScaleVec(const NOX::LAPACK::Vector& s) {
  scaleVec = s;
}

const NOX::Abstract::Vector&
LOCA::LAPACK::Group::getScaleVec() const {
  return scaleVec;
}

void
LOCA::LAPACK::Group::computeScaleVec() {
  scaleVec.init(1.0);
}

NOX::Abstract::Group::ReturnType 
LOCA::LAPACK::Group::augmentJacobianForHomotopy(double conParamValue)
{
  int size = scaleVec.length();

  // Scale the matrix by the value of the homotopy continuation param
  jacobianMatrix.scale(conParamValue);

  // Add the scaled identity matrix to the jacobian
  for (int i = 0; i < size; i++) 
    jacobianMatrix(i,i) += (1.0 - conParamValue);

  return NOX::Abstract::Group::Ok;
}
