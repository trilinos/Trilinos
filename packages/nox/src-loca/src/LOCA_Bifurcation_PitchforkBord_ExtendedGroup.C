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

#include "LOCA_Bifurcation_PitchforkBord_ExtendedGroup.H"
#include "LOCA_Bifurcation_TPBord_AbstractGroup.H"
#include "LOCA_Parameter_Vector.H"
#include "NOX_Parameter_List.H"

LOCA::Bifurcation::PitchforkBord::ExtendedGroup::ExtendedGroup(
                         LOCA::Bifurcation::TPBord::AbstractGroup& g,
			 const NOX::Abstract::Vector& asymVec,
			 const NOX::Abstract::Vector& lenVec,
			 int paramId)
  : grpPtr(&g), 
    pfXVec(g.getX(), lenVec, 0.0, 0.0),
    pfFVec(lenVec, lenVec, 0.0, 0.0),
    pfNewtonVec(lenVec, lenVec, 0.0, 0.0),
    asymVecPtr(asymVec.clone(NOX::DeepCopy)),
    lengthVecPtr(lenVec.clone(NOX::DeepCopy)), 
    bifParamId(paramId), 
    derivResidualParamPtr(lenVec.clone(NOX::ShapeCopy)), 
    derivNullResidualParamPtr(lenVec.clone(NOX::ShapeCopy)), 
    ownsGroup(false),
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false)
{
  init();
}

LOCA::Bifurcation::PitchforkBord::ExtendedGroup::ExtendedGroup(
                         const LOCA::Bifurcation::TPBord::AbstractGroup& g,
			 const NOX::Abstract::Vector& asymVec,
			 const NOX::Abstract::Vector& lenVec,
			 int paramId)
  : grpPtr(dynamic_cast<LOCA::Bifurcation::TPBord::AbstractGroup*>(g.clone())),
    pfXVec(g.getX(), lenVec, 0.0, 0.0),
    pfFVec(lenVec, lenVec, 0.0, 0.0),
    pfNewtonVec(lenVec, lenVec, 0.0, 0.0),
    asymVecPtr(asymVec.clone(NOX::DeepCopy)),
    lengthVecPtr(lenVec.clone(NOX::DeepCopy)), 
    bifParamId(paramId), 
    derivResidualParamPtr(lenVec.clone(NOX::ShapeCopy)), 
    derivNullResidualParamPtr(lenVec.clone(NOX::ShapeCopy)), 
    ownsGroup(true),
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false)
{
  init();
}

LOCA::Bifurcation::PitchforkBord::ExtendedGroup::ExtendedGroup(
               const LOCA::Bifurcation::PitchforkBord::ExtendedGroup& source, 
	       NOX::CopyType type)
  :  grpPtr(dynamic_cast<LOCA::Bifurcation::TPBord::AbstractGroup*>(source.grpPtr->clone())),
    pfXVec(source.pfXVec, type),
    pfFVec(source.pfFVec, type),
    pfNewtonVec(source.pfNewtonVec, type),
    asymVecPtr(source.asymVecPtr->clone(type)),
    lengthVecPtr(source.lengthVecPtr->clone(type)), 
    bifParamId(source.bifParamId),
    derivResidualParamPtr(source.derivResidualParamPtr->clone(type)),
    derivNullResidualParamPtr(source.derivNullResidualParamPtr->clone(type)),
    ownsGroup(true),
    isValidF(source.isValidF),
    isValidJacobian(source.isValidJacobian),
    isValidNewton(source.isValidNewton) {}


LOCA::Bifurcation::PitchforkBord::ExtendedGroup::~ExtendedGroup() 
{
  if (ownsGroup)
    delete grpPtr;
  delete asymVecPtr;
  delete lengthVecPtr;
  delete derivResidualParamPtr;
  delete derivNullResidualParamPtr;
}

LOCA::Continuation::AbstractGroup&
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::operator=(
                             const LOCA::Continuation::AbstractGroup& source)
{
  return *this = 
    dynamic_cast<const LOCA::Bifurcation::PitchforkBord::ExtendedGroup&>(source);
}

NOX::Abstract::Group&
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::operator=(
                                           const NOX::Abstract::Group& source)
{
  return *this = 
    dynamic_cast<const LOCA::Bifurcation::PitchforkBord::ExtendedGroup&>(source);
}

LOCA::Bifurcation::PitchforkBord::ExtendedGroup&
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::operator=(
                       const LOCA::Bifurcation::PitchforkBord::ExtendedGroup& source) 
{

  // Protect against A = A
  if (this != &source) {
    NOX::CopyType type = NOX::DeepCopy;

    // Delete old values
    delete asymVecPtr;
    delete lengthVecPtr;
    delete derivResidualParamPtr;
    delete derivNullResidualParamPtr;

    // Copy values
    *grpPtr = *source.grpPtr;
    pfXVec = source.pfXVec;
    pfFVec = source.pfFVec;
    pfNewtonVec = source.pfNewtonVec;
    asymVecPtr = source.asymVecPtr->clone(type);
    lengthVecPtr = source.lengthVecPtr->clone(type);
    derivResidualParamPtr = source.derivResidualParamPtr->clone(type);
    derivNullResidualParamPtr = source.derivNullResidualParamPtr->clone(type);
    bifParamId = source.bifParamId;
    isValidF = source.isValidF;
    isValidJacobian = source.isValidJacobian;
    isValidNewton = source.isValidNewton;
  }

  return *this;
}

NOX::Abstract::Group*
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::clone(NOX::CopyType type) const 
{
  return new LOCA::Bifurcation::PitchforkBord::ExtendedGroup(*this, type);
}

void
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::setParams(
                                               const LOCA::ParameterVector& p) 
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  grpPtr->setParams(p);
}

const LOCA::ParameterVector&
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::getParams() const 
{
  return grpPtr->getParams();
}

void
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::setParam(int paramID, double val)
{
  grpPtr->setParam(paramID, val);
}

double
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::getParam(int paramID) const
{
  return grpPtr->getParam(paramID);
}

void
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::setParam(string paramID, double val)
{
  grpPtr->setParam(paramID, val);
}

double
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::getParam(string paramID) const
{
  return grpPtr->getParam(paramID);
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::computeDfDp(int paramID, 
					      NOX::Abstract::Vector& result)
{
  NOX::Abstract::Group::ReturnType res;

  // Cast result to TP vector
  LOCA::Bifurcation::PitchforkBord::ExtendedVector& pf_result = 
    dynamic_cast<LOCA::Bifurcation::PitchforkBord::ExtendedVector&>(result);

  // Compute f, J, J*n
  res = computeF();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute df/dp
  res = grpPtr->computeDfDp(paramID, pf_result.getXVec());
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute d(Jn)/dp
  res = grpPtr->computeDJnDp(pfXVec.getNullVec(), paramID,
			     pfFVec.getNullVec(), pf_result.getNullVec());
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Set slack componenet
  pf_result.getSlackVar() = 0.0;

  // Set parameter componenet
  pf_result.getBifParam() = 0.0;

  return res;
}

double
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::getBifParam() const 
{
  return grpPtr->getParam(bifParamId);
}

void
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::setBifParam(double param) 
{
  grpPtr->setParam(bifParamId, param);

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

void
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::setX(const NOX::Abstract::Vector& y) 
{
  setX( dynamic_cast<const LOCA::Bifurcation::PitchforkBord::ExtendedVector&>(y) );
}

void
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::setX(
                   const LOCA::Bifurcation::PitchforkBord::ExtendedVector& y) 
{
  grpPtr->setX( y.getXVec() );
  pfXVec = y;
  setBifParam(pfXVec.getBifParam());

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

void
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::computeX(
					      const NOX::Abstract::Group& g, 
					      const NOX::Abstract::Vector& d,
					      double step) 
{
  computeX( dynamic_cast<const LOCA::Bifurcation::PitchforkBord::ExtendedGroup&>(g),
	    dynamic_cast<const LOCA::Bifurcation::PitchforkBord::ExtendedVector&>(d),
	    step);
}

void
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::computeX(
                            const LOCA::Bifurcation::PitchforkBord::ExtendedGroup& g, 
			    const LOCA::Bifurcation::PitchforkBord::ExtendedVector& d,
			    double step) 
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  grpPtr->computeX(*(g.grpPtr), d.getXVec(), step);
  pfXVec.update(1.0, g.getX(), step, d, 0.0);
  setBifParam(pfXVec.getBifParam());
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::computeF() 
{
  if (isValidF)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType res;

  const NOX::Abstract::Vector& x_x = pfXVec.getXVec();
  const NOX::Abstract::Vector& x_null = pfXVec.getNullVec();
  double x_slack = pfXVec.getSlackVar();

  NOX::Abstract::Vector& f_x = pfFVec.getXVec();
  NOX::Abstract::Vector& f_null = pfFVec.getNullVec();
  double& f_slack = pfFVec.getSlackVar();
  double& f_param = pfFVec.getBifParam();

  res = grpPtr->computeF();
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  f_x.update(1.0, grpPtr->getF(), x_slack, *asymVecPtr, 0.0);
  
  res = grpPtr->computeJacobian();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  res = grpPtr->applyJacobian(x_null, f_null);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  f_slack = grpPtr->innerProduct(x_x, *asymVecPtr);
  
  f_param = x_null.dot(*lengthVecPtr) - 1.0;
  
  isValidF = true;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::computeJacobian() 
{
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType res = computeF();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  res = grpPtr->computeJacobian();
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  res = grpPtr->computeDfDp(bifParamId, *derivResidualParamPtr);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  res = grpPtr->computeDJnDp(pfXVec.getNullVec(), bifParamId,
			     pfFVec.getNullVec(), *derivNullResidualParamPtr);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  isValidJacobian = true;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::computeGradient() 
{
  return NOX::Abstract::Group::NotDefined;
}
   
NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::computeNewton(
                                                 NOX::Parameter::List& params) 
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType res = computeF();
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  res = computeJacobian();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  res = applyJacobianInverse(params, pfFVec, pfNewtonVec);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  pfNewtonVec.scale(-1.0);
  isValidNewton = true;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::applyJacobian(
                                          const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  // Cast vectors to PitchforkBord::Vectors
  const LOCA::Bifurcation::PitchforkBord::ExtendedVector& pf_input = 
    dynamic_cast<const LOCA::Bifurcation::PitchforkBord::ExtendedVector&>(input);
  LOCA::Bifurcation::PitchforkBord::ExtendedVector& pf_result = 
    dynamic_cast<LOCA::Bifurcation::PitchforkBord::ExtendedVector&>(result);

  // Get constant references to input vector components
  const NOX::Abstract::Vector& input_x = pf_input.getXVec();
  const NOX::Abstract::Vector& input_null = pf_input.getNullVec();
  double input_slack = pf_input.getSlackVar();
  double input_param = pf_input.getBifParam();

  // Get non-constant references to result vector components
  NOX::Abstract::Vector& result_x = pf_result.getXVec();
  NOX::Abstract::Vector& result_null = pf_result.getNullVec();
  double& result_slack = pf_result.getSlackVar();
  double& result_param = pf_result.getBifParam();

  // Temporary vector
  NOX::Abstract::Vector *tmp = input_null.clone(NOX::ShapeCopy);

  // Return value
  NOX::Abstract::Group::ReturnType res;

  // compute J*x
  res = grpPtr->applyJacobian(input_x, result_x);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // compute J*x + sigma*psi + p*dR/dp
  result_x.update(input_param, *derivResidualParamPtr, input_slack, 
		  *asymVecPtr, 1.0);

  // compute J*y
  res = grpPtr->applyJacobian(input_null, result_null);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // compute J*y + p*dJy/dp
  result_null.update(input_param, *derivNullResidualParamPtr, 1.0);

  // compute (dJy/dx)*x
  res = grpPtr->computeDJnDxa(pfXVec.getNullVec(), input_x, 
			      pfFVec.getNullVec(), *tmp);
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  // compute (dJy/dx)*x + J*y + p*dJy/dp
  result_null.update(1.0, *tmp, 1.0);

  // compute <x,psi>
  result_slack = grpPtr->innerProduct(input_x, *asymVecPtr);

  // compute phi^T*y
  result_param = lengthVecPtr->dot(input_null);

  delete tmp;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::applyJacobianTranspose(
					  const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::applyJacobianInverse(
					 NOX::Parameter::List& params,
					 const NOX::Abstract::Vector& input,
					 NOX::Abstract::Vector& result) const 
{
  const NOX::Abstract::Vector* inputs[1];
  NOX::Abstract::Vector* results[1];

  inputs[0] = &input;
  results[0] = &result;

  return applyJacobianInverseMulti(params, inputs, results, 1);
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::applyJacobianInverseMulti(
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
    new const NOX::Abstract::Vector*[m+2];
  const NOX::Abstract::Vector** inputs_null =
    new const NOX::Abstract::Vector*[m+2];
  double* inputs_slacks = new double[m];
  double* inputs_params = new double[m];

  NOX::Abstract::Vector** results_x = new NOX::Abstract::Vector*[m+2];
  NOX::Abstract::Vector** results_null = new NOX::Abstract::Vector*[m+2];
  double** results_slacks = new double*[m];
  double** results_params = new double*[m];
  NOX::Abstract::Vector** tmp = new NOX::Abstract::Vector*[m+2];

  const LOCA::Bifurcation::PitchforkBord::ExtendedVector* constPitchforkVecPtr;
  LOCA::Bifurcation::PitchforkBord::ExtendedVector* pitchforkVecPtr;

  for (int i=0; i<m; i++) {
    constPitchforkVecPtr = 
      dynamic_cast<const LOCA::Bifurcation::PitchforkBord::ExtendedVector*>(inputs[i]);
    inputs_x[i] = &(constPitchforkVecPtr->getXVec());
    inputs_null[i] = &(constPitchforkVecPtr->getNullVec());
    inputs_slacks[i] = constPitchforkVecPtr->getSlackVar();
    inputs_params[i] = constPitchforkVecPtr->getBifParam();
 
    pitchforkVecPtr = 
      dynamic_cast<LOCA::Bifurcation::PitchforkBord::ExtendedVector*>(results[i]);
    results_x[i] = &(pitchforkVecPtr->getXVec());
    results_null[i] = &(pitchforkVecPtr->getNullVec());
    results_slacks[i] = &(pitchforkVecPtr->getSlackVar());
    results_params[i] = &(pitchforkVecPtr->getBifParam());

    tmp[i] = inputs_x[i]->clone(NOX::ShapeCopy);
  }

  // Set next to last components to deriv. w.r.t. parameter
  inputs_x[m] = derivResidualParamPtr;
  inputs_null[m] = derivNullResidualParamPtr;
  results_x[m] = inputs_x[m]->clone(NOX::ShapeCopy);
  results_null[m] = inputs_x[m]->clone(NOX::ShapeCopy);
  tmp[m] = inputs_x[m]->clone(NOX::ShapeCopy);

  inputs_x[m+1] = asymVecPtr;
  inputs_null[m+1] = NULL;  // really all zeros
  results_x[m+1] = inputs_x[m+1]->clone(NOX::ShapeCopy);
  results_null[m+1] = inputs_x[m+1]->clone(NOX::ShapeCopy);
  tmp[m+1] = inputs_x[m+1]->clone(NOX::ShapeCopy);

  // Solve J*results_x = inputs_x
  res = grpPtr->applySingularJacobianInverseMulti(params, inputs_x, 
						  pfXVec.getNullVec(),
						  pfFVec.getNullVec(),
						  results_x, m+2);
  if (res != NOX::Abstract::Group::Ok)
    return res;
  

  // Compute tmp = inputs_null - (dJy/dx)*results_x
  for (int i=0; i<m+2; i++) {
    res = grpPtr->computeDJnDxa(pfXVec.getNullVec(), *results_x[i],
				pfFVec.getNullVec(), *tmp[i]);

    if (res != NOX::Abstract::Group::Ok)
      return res;

    // inputs_null[m+1] is zero so don't do that update
    if (i == m+1)
      tmp[i]->scale(-1.0);
    else
      tmp[i]->update(1.0, *inputs_null[i], -1.0);
  } 

  // Solve J*results_null = tmp
  res = grpPtr->applySingularJacobianInverseMulti(params, tmp, 
						  pfXVec.getNullVec(),
						  pfFVec.getNullVec(),
						  results_null, m+2);
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  // Compute and set results
  double ipb = grpPtr->innerProduct(*results_x[m],   *asymVecPtr);
  double ipc = grpPtr->innerProduct(*results_x[m+1], *asymVecPtr);
  double lte = lengthVecPtr->dot(*results_null[m]);
  double ltf = lengthVecPtr->dot(*results_null[m+1]);
  double denom_slack = lte*ipc - ltf*ipb;
  double ipa, ltd;
 
  for (int i=0; i<m; i++) {
    ipa = grpPtr->innerProduct(*results_x[i], *asymVecPtr);
    ltd = lengthVecPtr->dot(*results_null[i]);
    *results_slacks[i] = (lte*(ipa - inputs_slacks[i]) - 
			 ipb*(ltd - inputs_params[i])) / denom_slack;
    *results_params[i] = (ltd - inputs_params[i] - ltf*(*results_slacks[i]))
      /lte;

    results_x[i]->update(-*results_params[i], *results_x[m], 
			 -*results_slacks[i], *results_x[m+1], 1.0);
    results_null[i]->update(-*results_params[i], *results_null[m], 
			    -*results_slacks[i], *results_null[m+1], 1.0);

    delete tmp[i];
  }

  delete results_x[m];
  delete results_null[m];
  delete tmp[m];

  delete results_x[m+1];
  delete results_null[m+1];
  delete tmp[m+1];

  delete [] inputs_x;
  delete [] inputs_null;
  delete [] inputs_params;
  delete [] results_x;
  delete [] results_null;
  delete [] results_slacks;
  delete [] results_params;
  delete [] tmp;

  return res;
}

bool
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::isF() const 
{
  return isValidF;
}

bool
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::isJacobian() const 
{
  return isValidJacobian;
}

bool
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::isGradient() const 
{
  return false;
}

bool
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::isNewton() const 
{
  return isValidNewton;
}
  
const NOX::Abstract::Vector&
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::getX() const 
{
  return pfXVec;
}

const NOX::Abstract::Vector&
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::getF() const 
{
  return pfFVec;
}

double
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::getNormF() const 
{
  return pfFVec.norm();
}

const NOX::Abstract::Vector&
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::getGradient() const 
{
  cout << "ERROR: LOCA::Bifurcation::PitchforkBord::ExtendedGroup::getGradient() "
       << " - not implemented" << endl;
  throw "LOCA Error";
  return getNewton();
}

const NOX::Abstract::Vector&
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::getNewton() const 
{
  return pfNewtonVec;
}

double
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::getNormNewtonSolveResidual() const 
{
  LOCA::Bifurcation::PitchforkBord::ExtendedVector residual = pfFVec;
  
  NOX::Abstract::Group::ReturnType res = applyJacobian(pfNewtonVec, residual);
  if (res != NOX::Abstract::Group::Ok) {
    cout << "ERROR: applyJacobian() in getNormNewtonSolveResidual "
	 << " returned not ok" << endl;
    throw "LOCA Error";
    return 0.0;
  }

  residual.update(1.0, pfFVec, 1.0);
  return residual.norm();
}

void
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::printSolution(const double conParam) const 
{
  cout << "LOCA::Bifurcation::PitchforkBord::ExtendedGroup::printSolution\n";
  cout << "\tSlack variable sigma = " << pfXVec.getSlackVar() << endl;
  cout << "\tPrinting Solution Vector for conParam = " << conParam << endl;
  grpPtr->printSolution(conParam);

  cout << "\tPrinting Null Vector for bif param = " << getBifParam() << endl;
  grpPtr->printSolution(pfXVec.getNullVec(), pfXVec.getBifParam());
}

const  LOCA::Continuation::AbstractGroup&
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::getUnderlyingGroup() const
{
  return grpPtr->getUnderlyingGroup();
}

 LOCA::Continuation::AbstractGroup&
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::getUnderlyingGroup()
{
  return grpPtr->getUnderlyingGroup();
}

void
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::init()
{
  pfXVec.getBifParam() = getBifParam();

  // Rescale length vector so that the normalization condition is met
  double lVecDotNullVec;
  lVecDotNullVec = lengthVecPtr->dot(pfXVec.getNullVec());
  if (lVecDotNullVec == 0.0) {
    cout << "ERROR: LOCA::Bifurcation::PitchforkBord::ExtendedGroup::init\n"
         << "     : length vector can not have Norm zero " << endl;

    throw "LOCA Error";
  }
  cout << "\tIn LOCA::Bifurcation::PitchforkBord::ExtendedGroup::init, scaling lenVec by:" 
       << (1.0 / lVecDotNullVec) << endl;
  lengthVecPtr->scale(1.0 / lVecDotNullVec);
}
