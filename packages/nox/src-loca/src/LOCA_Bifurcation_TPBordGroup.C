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

#include "LOCA_Bifurcation_TPBordGroup.H"
#include "LOCA_Parameter_Vector.H"
#include "NOX_Parameter_List.H"

LOCA::Bifurcation::TPBordGroup::TPBordGroup(const Abstract::Group& g,
					    const NOX::Abstract::Vector& lenVec,
					    int paramId,
					    const DerivUtils& d)
  : grpPtr(dynamic_cast<LOCA::Abstract::Group*>(g.clone(NOX::DeepCopy))), 
    tpXVec(g.getX(), lenVec, 0.0),
    tpFVec(lenVec, lenVec, 0.0),
    tpNewtonVec(lenVec, lenVec, 0.0),
    lengthVecPtr(lenVec.clone(NOX::ShapeCopy)), 
    tpScaleVec(g.getScaleVec(), g.getScaleVec(), 1.0),
    bifParamId(paramId), 
    derivResidualParamPtr(lenVec.clone(NOX::ShapeCopy)), 
    derivNullResidualParamPtr(lenVec.clone(NOX::ShapeCopy)), 
    derivPtr(d.clone(NOX::DeepCopy)), 
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false)
{
  tpXVec.getBifParam() = getBifParam();
}

LOCA::Bifurcation::TPBordGroup::TPBordGroup(const LOCA::Bifurcation::TPBordGroup& source, 
					    NOX::CopyType type)
  : grpPtr(dynamic_cast<LOCA::Abstract::Group*>(source.grpPtr->clone(type))), 
    tpXVec(source.tpXVec, type),
    tpFVec(source.tpFVec, type),
    tpNewtonVec(source.tpNewtonVec, type),
    lengthVecPtr(source.lengthVecPtr->clone(type)), 
    tpScaleVec(source.tpScaleVec, type),
    bifParamId(source.bifParamId),
    derivResidualParamPtr(source.derivResidualParamPtr->clone(type)),
    derivNullResidualParamPtr(source.derivNullResidualParamPtr->clone(type)),
    derivPtr(source.derivPtr->clone(type)),
    isValidF(source.isValidF),
    isValidJacobian(source.isValidJacobian),
    isValidNewton(source.isValidNewton) {}


LOCA::Bifurcation::TPBordGroup::~TPBordGroup() 
{
  delete grpPtr;
  delete lengthVecPtr;
  delete derivResidualParamPtr;
  delete derivNullResidualParamPtr;
  delete derivPtr;
  
}

LOCA::Abstract::Group&
LOCA::Bifurcation::TPBordGroup::operator=(const LOCA::Abstract::Group& source)
{
  return *this = dynamic_cast<const LOCA::Bifurcation::TPBordGroup&>(source);
}

NOX::Abstract::Group&
LOCA::Bifurcation::TPBordGroup::operator=(const NOX::Abstract::Group& source)
{
  return *this = dynamic_cast<const LOCA::Bifurcation::TPBordGroup&>(source);
}

LOCA::Bifurcation::TPBordGroup&
LOCA::Bifurcation::TPBordGroup::operator=(const LOCA::Bifurcation::TPBordGroup& source) 
{

  // Protect against A = A
  if (this != &source) {
    NOX::CopyType type = NOX::DeepCopy;

    // Delete old values
    delete grpPtr;
    delete lengthVecPtr;
    delete derivResidualParamPtr;
    delete derivNullResidualParamPtr;
    delete derivPtr;

    // Copy values
    grpPtr = dynamic_cast<LOCA::Abstract::Group*>(source.grpPtr->clone(type));
    tpXVec = source.tpXVec;
    tpFVec = source.tpFVec;
    tpNewtonVec = source.tpNewtonVec;
    lengthVecPtr = source.lengthVecPtr->clone(type);
    tpScaleVec = source.tpScaleVec;
    derivResidualParamPtr = source.derivResidualParamPtr->clone(type);
    derivNullResidualParamPtr = source.derivNullResidualParamPtr->clone(type);
    derivPtr = source.derivPtr->clone(type);
    bifParamId = source.bifParamId;
    isValidF = source.isValidF;
    isValidJacobian = source.isValidJacobian;
    isValidNewton = source.isValidNewton;
  }

  return *this;
}

NOX::Abstract::Group*
LOCA::Bifurcation::TPBordGroup::clone(NOX::CopyType type) const 
{
  return new LOCA::Bifurcation::TPBordGroup(*this, type);
}

void
LOCA::Bifurcation::TPBordGroup::setParams(const LOCA::ParameterVector& p) 
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  grpPtr->setParams(p);
}

const LOCA::ParameterVector&
LOCA::Bifurcation::TPBordGroup::getParams() const 
{
  return grpPtr->getParams();
}

void
LOCA::Bifurcation::TPBordGroup::setParam(int paramID, double val)
{
  grpPtr->setParam(paramID, val);
}

double
LOCA::Bifurcation::TPBordGroup::getParam(int paramID) const
{
  return grpPtr->getParam(paramID);
}

void
LOCA::Bifurcation::TPBordGroup::setParam(string paramID, double val)
{
  grpPtr->setParam(paramID, val);
}

double
LOCA::Bifurcation::TPBordGroup::getParam(string paramID) const
{
  return grpPtr->getParam(paramID);
}

double
LOCA::Bifurcation::TPBordGroup::getBifParam() const 
{
  LOCA::ParameterVector params(grpPtr->getParams());
  return params[bifParamId];
}

void
LOCA::Bifurcation::TPBordGroup::setBifParam(double param) 
{
  LOCA::ParameterVector params(grpPtr->getParams());

  params[bifParamId] = param;
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  grpPtr->setParams(params);
}

void
LOCA::Bifurcation::TPBordGroup::setX(const NOX::Abstract::Vector& y) 
{
  setX( dynamic_cast<const LOCA::Bifurcation::TPBordVector&>(y) );
}

void
LOCA::Bifurcation::TPBordGroup::setX(const LOCA::Bifurcation::TPBordVector& y) 
{
  grpPtr->setX( y.getXVec() );
  tpXVec = y;

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

void
LOCA::Bifurcation::TPBordGroup::computeX(const NOX::Abstract::Group& g, 
					 const NOX::Abstract::Vector& d,
					 double step) 
{
  computeX( dynamic_cast<const LOCA::Bifurcation::TPBordGroup&>(g),
	    dynamic_cast<const LOCA::Bifurcation::TPBordVector&>(d),
	    step);
}

void
LOCA::Bifurcation::TPBordGroup::computeX(const LOCA::Bifurcation::TPBordGroup& g, 
					 const LOCA::Bifurcation::TPBordVector& d,
					 double step) 
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  grpPtr->computeX(*(g.grpPtr), d.getXVec(), step);
  tpXVec.update(1.0, g.getX(), step, d, 0.0);
  setBifParam(tpXVec.getBifParam());
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBordGroup::computeF() 
{
  if (isValidF)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType res;

  res = grpPtr->computeF();
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  tpFVec.getXVec() = grpPtr->getF();
  
  res = grpPtr->computeJacobian();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  res = grpPtr->applyJacobian(tpXVec.getNullVec(), 
			      tpFVec.getNullVec());
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  tpFVec.getBifParam() = tpXVec.getNullVec().dot(*lengthVecPtr) - 1.0;
  
  isValidF = true;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBordGroup::computeJacobian() 
{
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType res = computeF();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  res = grpPtr->computeJacobian();
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  res = derivPtr->computeDfDp(*grpPtr, bifParamId, 
			      *derivResidualParamPtr);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  res = derivPtr->computeDJnDp(*grpPtr, tpXVec.getNullVec(), 
			       bifParamId,
			       tpFVec.getNullVec(), 
			       *derivNullResidualParamPtr);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  isValidJacobian = true;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBordGroup::computeGradient() 
{
  return NOX::Abstract::Group::NotDefined;
}
   
NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBordGroup::computeNewton(NOX::Parameter::List& params) 
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType res = computeF();
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  res = computeJacobian();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  res = applyJacobianInverse(params, tpFVec, tpNewtonVec);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  tpNewtonVec.scale(-1.0);
  isValidNewton = true;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBordGroup::applyJacobian(const NOX::Abstract::Vector& input,
					      NOX::Abstract::Vector& result) const 
{
  // Cast vectors to TPBordVectors
  const LOCA::Bifurcation::TPBordVector& tp_input = 
    dynamic_cast<const LOCA::Bifurcation::TPBordVector&>(input);
  LOCA::Bifurcation::TPBordVector& tp_result = 
    dynamic_cast<LOCA::Bifurcation::TPBordVector&>(result);

  // Get constant references to input vector components
  const NOX::Abstract::Vector& input_x = tp_input.getXVec();
  const NOX::Abstract::Vector& input_null = tp_input.getNullVec();
  double input_param = tp_input.getBifParam();

  // Get non-constant references to result vector components
  NOX::Abstract::Vector& result_x = tp_result.getXVec();
  NOX::Abstract::Vector& result_null = tp_result.getNullVec();
  double& result_param = tp_result.getBifParam();

  // Temporary vector
  NOX::Abstract::Vector *tmp = input_null.clone(NOX::ShapeCopy);

  // Value of bifurcation parameter
  double bifParam = getBifParam();

  // Return value
  NOX::Abstract::Group::ReturnType res;

  // compute J*x
  res = grpPtr->applyJacobian(input_x, result_x);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // compute J*x + p*dR/dp
  result_x.update(bifParam, *derivResidualParamPtr, 1.0);

  // compute J*y
  res = grpPtr->applyJacobian(input_null, result_null);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // compute J*y + p*dJy/dp
  result_null.update(bifParam, *derivNullResidualParamPtr, 1.0);

  // compute (dJy/dx)*x
  res = derivPtr->computeDJnDxa(*grpPtr, tpXVec.getNullVec(), input_x, 
				tpFVec.getNullVec(), *tmp);
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  // compute (dJy/dx)*x + J*y + p*dJy/dp
  result_null.update(1.0, *tmp, 1.0);

  // compute phi^T*y
  result_param = lengthVecPtr->dot(input_null);

  delete tmp;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBordGroup::applyJacobianTranspose(const NOX::Abstract::Vector& input,
						       NOX::Abstract::Vector& result) const 
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBordGroup::applyJacobianInverse(NOX::Parameter::List& params,
						     const NOX::Abstract::Vector& input,
						     NOX::Abstract::Vector& result) const 
{
  // Cast vectors to TPBordVectors
  const LOCA::Bifurcation::TPBordVector& tp_input = 
    dynamic_cast<const LOCA::Bifurcation::TPBordVector&>(input);
  LOCA::Bifurcation::TPBordVector& tp_result = 
    dynamic_cast<LOCA::Bifurcation::TPBordVector&>(result);

  // Get constant references to input vector components
  const NOX::Abstract::Vector& input_x = tp_input.getXVec();
  const NOX::Abstract::Vector& input_null = tp_input.getNullVec();
  double input_param = tp_input.getBifParam();

  // Get non-constant references to result vector components
  NOX::Abstract::Vector& result_x = tp_result.getXVec();
  NOX::Abstract::Vector& result_null = tp_result.getNullVec();
  double& result_param = tp_result.getBifParam();

  // Tmporary vectors
  NOX::Abstract::Vector *a = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector *b = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector *c = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector *d = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector *derivJa = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector *derivJb = input_x.clone(NOX::ShapeCopy);
  
  // Return type
  NOX::Abstract::Group::ReturnType res;

  // Solve J*a = input_x
  // Solve J*b = dR/dp
  // res = underlyingGroupApplyJacobianInverseManager(params, input_x,
//             *derivResidualParamPtr, tpXVec.getNullVec(), input_null, *a, *b);
  res = underlyingGroupApplyJacobianInverseManager(params, input_x,
            *derivResidualParamPtr, tpXVec.getNullVec(), tpFVec.getNullVec(),
	    *a, *b);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute (dJy/dx)*a
  res = derivPtr->computeDJnDxa(*grpPtr, tpXVec.getNullVec(), *a, 
				tpFVec.getNullVec(), *derivJa);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute input_null - (dJy/dx)*a
  derivJa->update(1.0, input_null, -1.0);

  // Compute (dJy/dx)*b
  res = derivPtr->computeDJnDxa(*grpPtr, tpXVec.getNullVec(), *b, 
				tpFVec.getNullVec(), *derivJb);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute (dJy/dx)*b - dJy/dp
  derivJb->update(-1.0, *derivNullResidualParamPtr, 1.0);

  // Solve J*c = input_null - (dJy/dx)*a
  // Solve J*d = (dJy/dx)*b - dJy/dp
  // res = underlyingGroupApplyJacobianInverseManager(params, *derivJa, *derivJb,
//                                       tpXVec.getNullVec(), input_null, *c, *d);
  res = underlyingGroupApplyJacobianInverseManager(params, *derivJa, *derivJb,
                                      tpXVec.getNullVec(), tpFVec.getNullVec(),
				      *c, *d);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute result_param = -(input_param - phi^T*c)/(phi^T*d)
  result_param = (input_param- lengthVecPtr->dot(*c))/(lengthVecPtr->dot(*d));

  // Compute result_x = a - result_param*b 
  result_x.update(1.0, *a, -result_param, *b, 0.0);

  // Compute result_null = c + result_param*d
  result_null.update(1.0, *c, result_param, *d, 0.0);

  // Clean up memory
  delete a;
  delete b;
  delete c;
  delete d;
  delete derivJa;
  delete derivJb;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBordGroup::applyJacobianInverseMulti(
			    NOX::Parameter::List& params,
			    const NOX::Abstract::Vector* const* inputs,
			    NOX::Abstract::Vector** results, int nVecs) const 
{
  // Number of input vectors
  int m = nVecs; 

  // Return type
   NOX::Abstract::Group::ReturnType res;
  
  // Build arrays of solution, null vector and parameter components
  const NOX::Abstract::Vector** inputs_x = 
    new const NOX::Abstract::Vector*[m+1];
  const NOX::Abstract::Vector** inputs_null =
    new const NOX::Abstract::Vector*[m+1];
  double* inputs_params = new double[m];

  NOX::Abstract::Vector** tmp1 = new NOX::Abstract::Vector*[m+1];
  NOX::Abstract::Vector** tmp2 = new NOX::Abstract::Vector*[m+1];
  NOX::Abstract::Vector** tmp3 = new NOX::Abstract::Vector*[m+1];

  const LOCA::Bifurcation::TPBordVector* constTPVecPtr;

  for (int i=0; i<m; i++) {
    constTPVecPtr = 
      dynamic_cast<const LOCA::Bifurcation::TPBordVector*>(inputs[i]);
    inputs_x[i] = &(constTPVecPtr->getXVec());
    inputs_null[i] = &(constTPVecPtr->getNullVec());
    inputs_params[i] = constTPVecPtr->getBifParam();

    tmp1[i] = inputs_x[i]->clone(NOX::ShapeCopy);
    tmp2[i] = inputs_x[i]->clone(NOX::ShapeCopy);
    tmp3[i] = inputs_x[i]->clone(NOX::ShapeCopy);
  }

  // Set last components to deriv. w.r.t. parameter
  inputs_x[m] = derivResidualParamPtr;
  inputs_null[m] = derivNullResidualParamPtr;
  tmp1[m] = inputs_x[m]->clone(NOX::ShapeCopy);
  tmp2[m] = inputs_x[m]->clone(NOX::ShapeCopy);
  tmp3[m] = inputs_x[m]->clone(NOX::ShapeCopy);

  // Solve J*tmp1 = inputs_x
  //res = grpPtr->applyJacobianInverseMulti(params, inputs_x, tmp1, m+1);
  res = underlyingGroupApplyJacobianInverseManagerMulti(params, inputs_x, 
							tpXVec.getNullVec(),
							tpFVec.getNullVec(),
							tmp1, m+1);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute tmp2 = (dJy/dx)*tmp1 - inputs_null
  for (int i=0; i<m+1; i++) {
    res = derivPtr->computeDJnDxa(*grpPtr, tpXVec.getNullVec(), *tmp1[i],
				  tpFVec.getNullVec(), *tmp2[i]);

    if (res != NOX::Abstract::Group::Ok)
      return res;

    tmp2[i]->update(-1.0, *inputs_null[i], 1.0);
  } 

  // Solve J*tmp3 = tmp2
  //res = grpPtr->applyJacobianInverseMulti(params, tmp2, tmp3, m+1);
  res = underlyingGroupApplyJacobianInverseManagerMulti(params, tmp2, 
							tpXVec.getNullVec(),
							tpFVec.getNullVec(),
							tmp3, m+1);
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  // Compute and set results
  double denom = lengthVecPtr->dot(*tmp3[m]);
  double w;
  LOCA::Bifurcation::TPBordVector* tpVecPtr;
  for (int i=0; i<m; i++) {
    tpVecPtr = dynamic_cast<LOCA::Bifurcation::TPBordVector*>(results[i]);

    w = (inputs_params[i] + lengthVecPtr->dot(*tmp3[i]))/denom;
    (tpVecPtr->getXVec()).update(1.0, *tmp1[i], -w, *tmp1[m], 0.0);
    (tpVecPtr->getNullVec()).update(-1.0, *tmp3[i], w, *tmp3[m], 0.0);
    tpVecPtr->getBifParam() = w;

    delete tmp1[i];
    delete tmp2[i];
    delete tmp3[i];
  }

  delete tmp1[m];
  delete tmp2[m];
  delete tmp3[m];

  delete [] tmp1;
  delete [] tmp2;
  delete [] tmp3;
  delete [] inputs_x;
  delete [] inputs_null;
  delete [] inputs_params;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBordGroup::underlyingGroupApplyJacobianInverseManager(
                       NOX::Parameter::List& params,
                       const NOX::Abstract::Vector& input1,
                       const NOX::Abstract::Vector& input2,
                       const NOX::Abstract::Vector& approxNullVec,
                       const NOX::Abstract::Vector& jacTimesApproxNullVec,
                       NOX::Abstract::Vector& result1,
                       NOX::Abstract::Vector& result2) const
{

  NOX::Abstract::Group::ReturnType res;
  string singularSolveOption = params.getParameter("Bifurcation Solve", "Default");

  if (singularSolveOption == "Nic") {

    // Solve of near-singular matrix with Nic method
    res = grpPtr->applyJacobianInverseNic(params, input1, approxNullVec,
                                          jacTimesApproxNullVec, result1);
    if (res != NOX::Abstract::Group::Ok)
      return res;

    res = grpPtr->applyJacobianInverseNic(params, input2, approxNullVec,
                                          jacTimesApproxNullVec, result2);

  }
  else if (singularSolveOption == "NicDay") {

    // Solve of near-singular matrix with Nic method
    res = grpPtr->applyJacobianInverseNicDay(params, input1, approxNullVec,
                                             jacTimesApproxNullVec, result1);
    if (res != NOX::Abstract::Group::Ok)
      return res;

    res = grpPtr->applyJacobianInverseNicDay(params, input2, approxNullVec,
                                             jacTimesApproxNullVec, result2);
  }
  else if (singularSolveOption == "Iterative Refinement") {

    // Solve of near-singular matrix with Nic method
    res = grpPtr->applyJacobianInverseItRef(params, input1, result1);
    if (res != NOX::Abstract::Group::Ok)
      return res;

    res = grpPtr->applyJacobianInverseItRef(params, input2, result2);

  }
  else {
    if (singularSolveOption != "Default") {
      cout << "WARNING: LOCA::Bifurcation::TPBordGroup::"
           << "underlyingGroupApplyJacobianInverseManager\n\t"
           << "Unknown Bifurcation Solve Option<" << singularSolveOption
           << ">:  Resetting to Default." << endl;

      params.setParameter("Bifurcation Solve", "Default");
    }
    
    // Default solve of matrix that is approaching singular, sequential RHS
    res = grpPtr->applyJacobianInverse(params, input1, result1);
    if (res != NOX::Abstract::Group::Ok)
      return res;

    res = grpPtr->applyJacobianInverse(params, input2, result2);
  }

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBordGroup::underlyingGroupApplyJacobianInverseManagerMulti(
                       NOX::Parameter::List& params,
                       const NOX::Abstract::Vector*const* inputs,
                       const NOX::Abstract::Vector& approxNullVec,
                       const NOX::Abstract::Vector& jacTimesApproxNullVec,
                       NOX::Abstract::Vector** results, int nVecs) const
{

  NOX::Abstract::Group::ReturnType res;
  string singularSolveOption = params.getParameter("Bifurcation Solve", "Default");

  if (singularSolveOption == "Nic") {

    // Solve of near-singular matrix with Nic method
    res = grpPtr->applyJacobianInverseNicMulti(params, inputs, approxNullVec,
					       jacTimesApproxNullVec, results,
					       nVecs);
  }
  else if (singularSolveOption == "NicDay") {

    // Solve of near-singular matrix with Nic method
    res = grpPtr->applyJacobianInverseNicDayMulti(params, inputs, 
						  approxNullVec,
						  jacTimesApproxNullVec, 
						  results, nVecs);
  }
  else if (singularSolveOption == "Iterative Refinement") {

    // Solve of near-singular matrix with Nic method
    res = grpPtr->applyJacobianInverseItRefMulti(params, inputs, results, 
						 nVecs);
  }
  else {
    if (singularSolveOption != "Default") {
      cout << "WARNING: LOCA::Bifurcation::TPBordGroup::"
           << "underlyingGroupApplyJacobianInverseManager\n\t"
           << "Unknown Bifurcation Solve Option<" << singularSolveOption
           << ">:  Resetting to Default." << endl;

      params.setParameter("Bifurcation Solve", "Default");
    }
    
    // Default solve of matrix that is approaching singular
    res = grpPtr->applyJacobianInverseMulti(params, inputs, results, nVecs);
  }

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBordGroup::applyRightPreconditioning(NOX::Parameter::List& params, const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const 
{
  // Cast vectors to TPBordVectors
  const LOCA::Bifurcation::TPBordVector& tp_input = 
    dynamic_cast<const LOCA::Bifurcation::TPBordVector&>(input);
  LOCA::Bifurcation::TPBordVector& tp_result = 
    dynamic_cast<LOCA::Bifurcation::TPBordVector&>(result);

  // Get constant references to input vector components
  const NOX::Abstract::Vector& input_x = tp_input.getXVec();
  const NOX::Abstract::Vector& input_null = tp_input.getNullVec();
  double input_param = tp_input.getBifParam();

  // Get non-constant references to result vector components
  NOX::Abstract::Vector& result_x = tp_result.getXVec();
  NOX::Abstract::Vector& result_null = tp_result.getNullVec();
  double& result_param = tp_result.getBifParam();

  // Temporary vectors
  NOX::Abstract::Vector *a = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector *b = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector *c = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector *d = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector *derivJ = input_x.clone(NOX::ShapeCopy);

  // Return value
  NOX::Abstract::Group::ReturnType res;

  // Solve J*a = input_x
  res = grpPtr->applyRightPreconditioning(params, input_x, *a);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Solve J*b = dR/dp
  res = grpPtr->applyRightPreconditioning(params, *derivResidualParamPtr, *b);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute (dJy/dx)*a
  res = derivPtr->computeDJnDxa(*grpPtr, tpXVec.getNullVec(), *a, 
				tpFVec.getNullVec(), *derivJ);

  // Compute input_null - (dJy/dx)*a
  *derivJ = derivJ->update(1.0, input_null, -1.0);

  // Solve J*c = input_null - (dJy/dx)*a
  res = grpPtr->applyRightPreconditioning(params, *derivJ, *c);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute (dJy/dx)*b
  res = derivPtr->computeDJnDxa(*grpPtr, tpXVec.getNullVec(), *b, 
				tpFVec.getNullVec(), *derivJ);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute (dJy/dx)*b - dJy/dp
  *derivJ = derivJ->update(-1.0, *derivNullResidualParamPtr, 1.0);

  // Solve J*d = (dJy/dx)*b - dJy/dp
  res = grpPtr->applyRightPreconditioning(params, *derivJ, *d);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute result_param = -(input_param - phi^T*c)/(phi^T*d)
  result_param = (input_param- lengthVecPtr->dot(*c))/(lengthVecPtr->dot(*d));

  // Compute result_x = a - result_param*b 
  result_x.update(1.0, *a, -result_param, *b, 0.0);

  // Compute result_null = c + result_param*d
  result_null.update(1.0, *c, result_param, *d, 0.0);

  // Clean up memory
  delete a;
  delete b;
  delete c;
  delete d;
  delete derivJ;

  return res;
}

bool
LOCA::Bifurcation::TPBordGroup::isF() const 
{
  return isValidF;
}

bool
LOCA::Bifurcation::TPBordGroup::isJacobian() const 
{
  return isValidJacobian;
}

bool
LOCA::Bifurcation::TPBordGroup::isGradient() const 
{
  return false;
}

bool
LOCA::Bifurcation::TPBordGroup::isNewton() const 
{
  return isValidNewton;
}
  
const NOX::Abstract::Vector&
LOCA::Bifurcation::TPBordGroup::getX() const 
{
  return tpXVec;
}

const NOX::Abstract::Vector&
LOCA::Bifurcation::TPBordGroup::getF() const 
{
  return tpFVec;
}

double
LOCA::Bifurcation::TPBordGroup::getNormF() const 
{
  return tpFVec.norm();
}

const NOX::Abstract::Vector&
LOCA::Bifurcation::TPBordGroup::getGradient() const 
{
  cout << "ERROR: LOCA::Bifurcation::TPBordGroup::getGradient() "
       << " - not implemented" << endl;
  throw "LOCA Error";
  return getNewton();
}

const NOX::Abstract::Vector&
LOCA::Bifurcation::TPBordGroup::getNewton() const 
{
  return tpNewtonVec;
}

double
LOCA::Bifurcation::TPBordGroup::getNormNewtonSolveResidual() const 
{
  LOCA::Bifurcation::TPBordVector residual = tpFVec;
  
  NOX::Abstract::Group::ReturnType res = applyJacobian(tpNewtonVec, residual);
  if (res != NOX::Abstract::Group::Ok) {
    cout << "ERROR: applyJacobian() in getNormNewtonSolveResidual "
	 << " returned not ok" << endl;
    throw "LOCA Error";
    return 0.0;
  }

  residual.update(1.0, tpFVec, 1.0);
  return residual.norm();
}

void
LOCA::Bifurcation::TPBordGroup::print() const 
{
  cout << "Beginning TPBordGroup.print:" << endl;
  cout << "Group = " << endl;
  grpPtr->print();

  cout << endl;

  cout << "tpXVec = " << endl;
  tpXVec.print();

  cout << endl;

  if (isValidF) {
    cout << "tpFVec = " << endl;
    tpFVec.print();
  }
  else
    cout << "tpFVec not computed" << endl;

  cout << endl;

  if (isValidNewton) {
    cout << "tpNewtonVec = " << endl;
    tpNewtonVec.print();
  }
  else
    cout << "tpNewtonVec not computed" << endl;

  cout << endl;

  cout << "lengthVec = " << endl;
  lengthVecPtr->print();

  cout << endl;

  cout << "bifParmId = " << bifParamId << endl;

  cout << endl;

  if (isValidJacobian) {
    cout << "derivResdiualParam = " << endl;
    derivResidualParamPtr->print();
    cout << endl;
    cout << "derivNullResidualParam = " << endl;
    derivNullResidualParamPtr->print();
  }
  else
    cout << "deriv vectors not computed" << endl;

  cout << "End TPBordGroup.print:" << endl;
  cout << endl;
}

void
LOCA::Bifurcation::TPBordGroup::setScaleVec(const NOX::Abstract::Vector& s) {
  setScaleVec( dynamic_cast<const LOCA::Bifurcation::TPBordVector&>(s) );
}

void
LOCA::Bifurcation::TPBordGroup::setScaleVec(const LOCA::Bifurcation::TPBordVector& s) {
  tpScaleVec = s;
}

const NOX::Abstract::Vector&
LOCA::Bifurcation::TPBordGroup::getScaleVec() const {
  return tpScaleVec;
}
