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

LOCA::Abstract::Group::Group(const LOCA::DerivUtils& d)
  : derivPtr(d.clone(NOX::DeepCopy))
{
}

LOCA::Abstract::Group::Group(const LOCA::Abstract::Group& source, 
			     NOX::CopyType type)
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

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyJacobianInverseNic(
                          NOX::Parameter::List& params,
                          const NOX::Abstract::Vector& input,
                          const NOX::Abstract::Vector& approxNullVec,
                          const NOX::Abstract::Vector& JacApproxNullVec,
                          NOX::Abstract::Vector& result) const
{
  double alpha = approxNullVec.dot(input)
               / approxNullVec.dot(JacApproxNullVec);

  NOX::Abstract::Vector* tmpInput  = input.clone(NOX::DeepCopy);
  tmpInput->update(-alpha, JacApproxNullVec, 1.0);

  NOX::Abstract::Group::ReturnType
    res = applyJacobianInverse(params, *tmpInput, result);

  delete tmpInput;

  result.update(alpha, approxNullVec, 1.0);

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyJacobianInverseNicMulti(
                          NOX::Parameter::List& params,
                          const NOX::Abstract::Vector*const* inputs,
                          const NOX::Abstract::Vector& approxNullVec,
                          const NOX::Abstract::Vector& JacApproxNullVec,
                          NOX::Abstract::Vector** results, int nVecs) const
{
  double denom = approxNullVec.dot(JacApproxNullVec);

  double* alphas = new double[nVecs];
  NOX::Abstract::Vector** tmpInputs  = new NOX::Abstract::Vector*[nVecs];

  for (int i=0; i<nVecs; i++) {
    alphas[i] = approxNullVec.dot(*(inputs[i])) / denom;
    tmpInputs[i] = inputs[i]->clone(NOX::DeepCopy);
    tmpInputs[i]->update(-alphas[i], JacApproxNullVec, 1.0);
  }

  NOX::Abstract::Group::ReturnType
    res = applyJacobianInverseMulti(params, tmpInputs, results, nVecs);

  for (int i=0; i<nVecs; i++) {
    results[i]->update(alphas[i], approxNullVec, 1.0);
    delete tmpInputs[i];
  }

  delete [] tmpInputs;
  delete [] alphas;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyJacobianInverseNicDay(
                          NOX::Parameter::List& params,
                          const NOX::Abstract::Vector& input,
                          const NOX::Abstract::Vector& approxNullVec,
                          const NOX::Abstract::Vector& JacApproxNullVec,
                          NOX::Abstract::Vector& result) const
{
  double alpha = JacApproxNullVec.dot(input)
               / JacApproxNullVec.dot(JacApproxNullVec);

  NOX::Abstract::Vector* tmpInput  = input.clone(NOX::DeepCopy);
  tmpInput->update(-alpha, JacApproxNullVec, 1.0);

  NOX::Abstract::Group::ReturnType
    res = applyJacobianInverse(params, *tmpInput, result);

  delete tmpInput;

  result.update(alpha, approxNullVec, 1.0);

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyJacobianInverseNicDayMulti(
                          NOX::Parameter::List& params,
                          const NOX::Abstract::Vector*const* inputs,
                          const NOX::Abstract::Vector& approxNullVec,
                          const NOX::Abstract::Vector& JacApproxNullVec,
                          NOX::Abstract::Vector** results, int nVecs) const
{
  double denom = JacApproxNullVec.dot(JacApproxNullVec);

  double* alphas = new double[nVecs];
  NOX::Abstract::Vector** tmpInputs  = new NOX::Abstract::Vector*[nVecs];

  for (int i=0; i<nVecs; i++) {
    alphas[i] = JacApproxNullVec.dot(*(inputs[i]));
    tmpInputs[i] = inputs[i]->clone(NOX::DeepCopy);
    tmpInputs[i]->update(-alphas[i], JacApproxNullVec, 1.0);
  }

  NOX::Abstract::Group::ReturnType
    res = applyJacobianInverseMulti(params, tmpInputs, results, nVecs);

  for (int i=0; i<nVecs; i++) {
    results[i]->update(alphas[i], approxNullVec, 1.0);
    delete tmpInputs[i];
  }

  delete [] tmpInputs;
  delete [] alphas;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyJacobianInverseItRef(
                            NOX::Parameter::List& params,
                            const NOX::Abstract::Vector& input,
                            NOX::Abstract::Vector& result) const
{
  NOX::Abstract::Group::ReturnType 
    res = applyJacobianInverse(params, input, result);

  NOX::Abstract::Vector* remainder = input.clone(NOX::ShapeCopy);

  res = applyJacobian(result, *remainder);

  // r = b-Ax
  remainder->update(1.0, input, -1.0);

  NOX::Abstract::Vector* refinement = input.clone(NOX::ShapeCopy);

  // Ay=r
  res = applyJacobianInverse(params, *remainder, *refinement);

  // x+=y
  result.update(1.0, *refinement, 1.0);
  
  delete remainder;
  delete refinement;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyJacobianInverseItRefMulti(
                            NOX::Parameter::List& params,
                            const NOX::Abstract::Vector*const* inputs,
                            NOX::Abstract::Vector** results,
			    int nVecs) const
{
  NOX::Abstract::Vector** remainders = new NOX::Abstract::Vector*[nVecs];
  NOX::Abstract::Vector** refinements = new NOX::Abstract::Vector*[nVecs];

  NOX::Abstract::Group::ReturnType 
    res = applyJacobianInverseMulti(params, inputs, results, nVecs);

  for (int i=0; i<nVecs; i++) {
    remainders[i] = inputs[i]->clone(NOX::ShapeCopy);
    refinements[i] = inputs[i]->clone(NOX::ShapeCopy);

    res = applyJacobian(*(results[i]), *(remainders[i]));

    // r = b-Ax
    remainders[i]->update(1.0, *(inputs[i]), -1.0);
  }

  // Ay=r
  res = applyJacobianInverseMulti(params, remainders, refinements, nVecs);

  // x+=y
  for (int i=0; i<nVecs; i++) {
    results[i]->update(1.0, *(refinements[i]), 1.0);
    delete remainders[i];
    delete refinements[i];
  }
  
  delete [] remainders;
  delete [] refinements;

  return res;
}


NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyJacobianInverseMulti(NOX::Parameter::List& params,
			    const NOX::Abstract::Vector* const* inputs,
			    NOX::Abstract::Vector** outputs, int nVecs) const
{
  NOX::Abstract::Group::ReturnType res;

  for (int i=0; i<nVecs; i++) {
    res = applyJacobianInverse(params, *inputs[i], *outputs[i]);
    if (res != NOX::Abstract::Group::Ok)
      return res;
  }

  return res;
}

double
LOCA::Abstract::Group::innerProduct(const NOX::Abstract::Vector& x,
				    const NOX::Abstract::Vector& y) 
{
  return x.dot(y);
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::computeDfDp(int paramID, NOX::Abstract::Vector& result) 
{
  return derivPtr->computeDfDp(*this,paramID, result);
}
