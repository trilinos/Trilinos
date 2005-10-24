// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "LOCA_Bifurcation_PitchforkBord_ExtendedGroup.H"
#include "LOCA_Bifurcation_TPBord_AbstractGroup.H"
#include "LOCA_Parameter_Vector.H"
#include "NOX_Parameter_List.H"
#include "LOCA_Utils.H"
#include "LOCA_ErrorCheck.H"

LOCA::Bifurcation::PitchforkBord::ExtendedGroup::ExtendedGroup(
			      LOCA::Bifurcation::TPBord::AbstractGroup& g,
			      NOX::Parameter::List& bifParamList)
  : grpPtr(&g),
    pfXVec(g.getX(), g.getX(), 0.0, 0.0),
    pfFVec(g.getX(), g.getX(), 0.0, 0.0),
    pfNewtonVec(g.getX(), g.getX(), 0.0, 0.0),
    asymVecPtr(NULL),
    lengthVecPtr(NULL), 
    bifParamId(0), 
    derivResidualParamPtr(NULL), 
    derivNullResidualParamPtr(NULL), 
    ownsGroup(false),
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false)
{
  const char *func = "LOCA::Bifurcation::PitchforkBord::ExtendedGroup()";

  if (!bifParamList.isParameter("Bifurcation Parameter")) {
    LOCA::ErrorCheck::throwError(func,
				 "\"Bifurcation Parameter\" name is not set!");
  }
  string bifParamName = bifParamList.getParameter("Bifurcation Parameter",
						  "None");
  const ParameterVector& p = grpPtr->getParams();
  bifParamId = p.getIndex(bifParamName);

  if (!bifParamList.isParameter("Asymmetric Vector")) {
    LOCA::ErrorCheck::throwError(func,
				"\"Asymmetric Vector\" is not set!");
  }
  const NOX::Abstract::Vector* asVecPtr = 
    bifParamList.getAnyPtrParameter<NOX::Abstract::Vector>("Asymmetric Vector");

  if (!bifParamList.isParameter("Length Normalization Vector")) {
    LOCA::ErrorCheck::throwError(func,
				"\"Length Normalization Vector\" is not set!");
  }
  const NOX::Abstract::Vector* lenVecPtr = 
    bifParamList.getAnyPtrParameter<NOX::Abstract::Vector>("Length Normalization Vector");

  if (!bifParamList.isParameter("Initial Null Vector")) {
    LOCA::ErrorCheck::throwError(func,
				 "\"Initial Null Vector\" is not set!");
  }
  const NOX::Abstract::Vector* nullVecPtr = 
    bifParamList.getAnyConstPtrParameter<NOX::Abstract::Vector>("Initial Null Vector");

  bool perturbSoln = bifParamList.getParameter("Perturb Initial Solution", 
					       false);
  double perturbSize = bifParamList.getParameter("Relative Perturbation Size", 
						 1.0e-3);

  asymVecPtr = asVecPtr->clone(NOX::DeepCopy);
  lengthVecPtr = lenVecPtr->clone(NOX::DeepCopy);
  derivResidualParamPtr = lenVecPtr->clone(NOX::ShapeCopy);
  derivNullResidualParamPtr = lenVecPtr->clone(NOX::ShapeCopy);
  pfXVec.getNullVec() = *nullVecPtr;

  init(perturbSoln, perturbSize);
}

LOCA::Bifurcation::PitchforkBord::ExtendedGroup::ExtendedGroup(
                         LOCA::Bifurcation::TPBord::AbstractGroup& g,
			 const NOX::Abstract::Vector& asymVec,
			 const NOX::Abstract::Vector& lenVec,
			 const NOX::Abstract::Vector& nullVec,
			 int paramId)
  : grpPtr(&g), 
    pfXVec(g.getX(), nullVec, 0.0, 0.0),
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
			 const NOX::Abstract::Vector& nullVec,
			 int paramId)
  : grpPtr(dynamic_cast<LOCA::Bifurcation::TPBord::AbstractGroup*>(g.clone())),
    pfXVec(g.getX(), nullVec, 0.0, 0.0),
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

LOCA::Extended::AbstractGroup&
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::operator=(
                             const LOCA::Extended::AbstractGroup& source)
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

    // Copy values
    *grpPtr = *source.grpPtr;
    pfXVec = source.pfXVec;
    pfFVec = source.pfFVec;
    pfNewtonVec = source.pfNewtonVec;
    *asymVecPtr = *source.asymVecPtr;
    *lengthVecPtr = *source.lengthVecPtr;
    *derivResidualParamPtr = *source.derivResidualParamPtr;
    *derivNullResidualParamPtr = *source.derivNullResidualParamPtr;
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
  string callingFunction = 
    "LOCA::Bifurcation::PitchforkBord::ExtendedGroup::computeDfDp()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Cast result to pitchfork vector
  LOCA::Bifurcation::PitchforkBord::ExtendedVector& pf_result = 
    dynamic_cast<LOCA::Bifurcation::PitchforkBord::ExtendedVector&>(result);

  // Make sure F is valid
  if (!isF()) {
    status = computeF();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Compute df/dp
  status = grpPtr->computeDfDp(paramID, pf_result.getXVec());
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute d(Jn)/dp
  status = grpPtr->computeDJnDp(pfXVec.getNullVec(), paramID,
			     pfFVec.getNullVec(), pf_result.getNullVec());
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Set slack componenet
  pf_result.getSlackVar() = 0.0;

  // Set parameter componenet
  pf_result.getBifParam() = 0.0;

  return finalStatus;
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
  
  string callingFunction = 
    "LOCA::Bifurcation::PitchforkBord::ExtendedGroup::computeF()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  const NOX::Abstract::Vector& x_x = pfXVec.getXVec();
  const NOX::Abstract::Vector& x_null = pfXVec.getNullVec();
  double x_slack = pfXVec.getSlackVar();

  NOX::Abstract::Vector& f_x = pfFVec.getXVec();
  NOX::Abstract::Vector& f_null = pfFVec.getNullVec();
  double& f_slack = pfFVec.getSlackVar();
  double& f_param = pfFVec.getBifParam();

  // Compute underlying F
  if (!grpPtr->isF()) {
    status = grpPtr->computeF();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }
  f_x.update(1.0, grpPtr->getF(), x_slack, *asymVecPtr, 0.0);

  // Compute underlying Jacobian
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }
  
  // Compute J*n
  status = grpPtr->applyJacobian(x_null, f_null);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute <x,psi>
  f_slack = grpPtr->innerProduct(x_x, *asymVecPtr);
  
  // compute phi^T*n
  f_param = lTransNorm(x_null) - 1.0;
  
  isValidF = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::computeJacobian() 
{
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::Bifurcation::PitchforkBord::ExtendedGroup::computeJacobian()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // We need data from computeF() here, so make sure it's valid
  if (!isF()) {
    status = computeF();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Compute underlying df/dp (may invalidate underlying data)
  status = grpPtr->computeDfDp(bifParamId, *derivResidualParamPtr);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute underlying dJn/dp
  status = grpPtr->computeDJnDp(pfXVec.getNullVec(), bifParamId,
				pfFVec.getNullVec(), 
				*derivNullResidualParamPtr);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute underlying Jacobian
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }
  
  isValidJacobian = true;

  return finalStatus;
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

  string callingFunction = 
    "LOCA::Bifurcation::PitchforkBord::ExtendedGroup::computeNewton()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Make sure F is valid
  if (!isF()) {
    status = computeF();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }
  
  // Make sure Jacobian is valid
  if (!isJacobian()) {
    status = computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // zero out newton vec -- used as initial guess for some linear solvers
  pfNewtonVec.init(0.0);

  status = applyJacobianInverse(params, pfFVec, pfNewtonVec);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  pfNewtonVec.scale(-1.0);
  isValidNewton = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::applyJacobian(
                                          const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  string callingFunction = 
    "LOCA::Bifurcation::PitchforkBord::ExtendedGroup::applyJacobian()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  if (!isJacobian()) {
    LOCA::ErrorCheck::throwError(callingFunction,
				 "Called with invalid Jacobian!");
  }

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

  // verify underlying Jacobian is valid
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // compute J*x
  status = grpPtr->applyJacobian(input_x, result_x);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // compute J*x + sigma*psi + p*dR/dp
  result_x.update(input_param, *derivResidualParamPtr, input_slack, 
		  *asymVecPtr, 1.0);

  // compute J*y
  status = grpPtr->applyJacobian(input_null, result_null);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // compute J*y + p*dJy/dp
  result_null.update(input_param, *derivNullResidualParamPtr, 1.0);

  // compute (dJy/dx)*x
  status = grpPtr->computeDJnDxa(pfXVec.getNullVec(), input_x, 
				 pfFVec.getNullVec(), *tmp);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  
  // compute (dJy/dx)*x + J*y + p*dJy/dp
  result_null.update(1.0, *tmp, 1.0);

  // compute <x,psi>
  result_slack = grpPtr->innerProduct(input_x, *asymVecPtr);

  // compute phi^T*y
  result_param = lTransNorm(input_null);

  delete tmp;

  return finalStatus;
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
  string callingFunction = 
"LOCA::Bifurcation::PitchforkBord::ExtendedGroup::applyJacobianInverseMulti()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  if (!isJacobian()) {
    LOCA::ErrorCheck::throwError(callingFunction,
				 "Called with invalid Jacobian!");
  }

  // Number of input vectors
  int m = nVecs; 
  
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

  // verify underlying Jacobian is valid
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Solve J*results_x = inputs_x
  status = grpPtr->applySingularJacobianInverseMulti(params, inputs_x, 
						     pfXVec.getNullVec(),
						     pfFVec.getNullVec(),
						     results_x, m+2);
  finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);

  // Compute tmp = inputs_null - (dJy/dx)*results_x
  for (int i=0; i<m+2; i++) {
    status = grpPtr->computeDJnDxa(pfXVec.getNullVec(), *results_x[i],
				   pfFVec.getNullVec(), *tmp[i]);

    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);

    // inputs_null[m+1] is zero so don't do that update
    if (i == m+1)
      tmp[i]->scale(-1.0);
    else
      tmp[i]->update(1.0, *inputs_null[i], -1.0);
  }

  // verify underlying Jacobian is valid
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Solve J*results_null = tmp
  status = grpPtr->applySingularJacobianInverseMulti(params, tmp, 
						     pfXVec.getNullVec(),
						     pfFVec.getNullVec(),
						     results_null, m+2);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  
  // Compute and set results
  double ipb = grpPtr->innerProduct(*results_x[m],   *asymVecPtr);
  double ipc = grpPtr->innerProduct(*results_x[m+1], *asymVecPtr);
  double lte = lTransNorm(*results_null[m]);
  double ltf = lTransNorm(*results_null[m+1]);
  double denom_slack = lte*ipc - ltf*ipb;
  double ipa, ltd;
 
  for (int i=0; i<m; i++) {
    ipa = grpPtr->innerProduct(*results_x[i], *asymVecPtr);
    ltd = lTransNorm(*results_null[i]);
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

  return finalStatus;
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
  LOCA::ErrorCheck::throwError(
	     "LOCA::Bifurcation::PitchforkBord::ExtendedGroup::getGradient()",
	     " - not implemented");
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
  string callingFunction = 
    "LOCA::Bifurcation::PitchforkBord::ExtendedGroup::getNormNewtonSolveResidual()";
  NOX::Abstract::Group::ReturnType finalStatus;
  LOCA::Bifurcation::PitchforkBord::ExtendedVector residual = pfFVec;
  
  finalStatus = applyJacobian(pfNewtonVec, residual);
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  residual.update(1.0, pfFVec, 1.0);
  return residual.norm();
}

void
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::printSolution(const double conParam) const 
{
  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << endl
	 << "LOCA::Bifurcation::PitchforkBord::ExtendedGroup::printSolution\n";
    cout << "\tSlack variable sigma = " 
	 << LOCA::Utils::sci(pfXVec.getSlackVar()) << endl;
    cout << "\tPrinting Solution Vector for conParam = " 
	 << LOCA::Utils::sci(conParam) << endl;
  }
  grpPtr->printSolution(conParam);

  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "\tPrinting Null Vector for bif param = " 
	 << LOCA::Utils::sci(getBifParam()) << endl;
  }
  grpPtr->printSolution(pfXVec.getNullVec(), pfXVec.getBifParam());
}

void
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::printSolution(
					   const NOX::Abstract::Vector& x_,
					   const double conParam) const 
{
  const LOCA::Bifurcation::PitchforkBord::ExtendedVector& pf_x = 
    dynamic_cast<const LOCA::Bifurcation::PitchforkBord::ExtendedVector&>(x_);

  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << endl
	 << "LOCA::Bifurcation::PitchforkBord::ExtendedGroup::printSolution\n";
    cout << "\tSlack variable sigma = " 
	 << LOCA::Utils::sci(pf_x.getSlackVar()) << endl;
    cout << "\tPrinting Solution Vector for conParam = " 
	 << LOCA::Utils::sci(conParam) << endl;
  }
  grpPtr->printSolution(pf_x.getXVec(), conParam);

  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "\tPrinting Null Vector for bif param = " 
	 << LOCA::Utils::sci(pf_x.getBifParam()) << endl;
  }
  grpPtr->printSolution(pf_x.getNullVec(), pf_x.getBifParam());
}

const  LOCA::Continuation::AbstractGroup&
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::getUnderlyingGroup() const
{
  return *grpPtr;
}

 LOCA::Continuation::AbstractGroup&
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::getUnderlyingGroup()
{
  return *grpPtr;
}

void
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::init(bool perturbSoln, 
						      double perturbSize)
{
  pfXVec.getBifParam() = getBifParam();

  // Rescale length vector so that the normalization condition is met
  double lVecDotNullVec;
  lVecDotNullVec = lTransNorm(pfXVec.getNullVec());

  if (lVecDotNullVec == 0.0) {
    LOCA::ErrorCheck::throwError(
		  "LOCA::Bifurcation::PitchforkBord::ExtendedGroup::init()",
		  "null vector can be orthogonal to length-scaling vector");
  }
  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "\tIn LOCA::Bifurcation::PitchforkBord::ExtendedGroup::init, scaling null vector by:" 
	 << LOCA::Utils::sci(1.0 / lVecDotNullVec) << endl;
  }
  pfXVec.getNullVec().scale(1.0 / lVecDotNullVec);

  // Rescale asymmetric vector to have unit length
  double psi_norm = sqrt( grpPtr->innerProduct(*asymVecPtr,*asymVecPtr) );
  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "\tIn LOCA::Bifurcation::PitchforkBord::ExtendedGroup::init, scaling asymmetric vector by:" 
	 << LOCA::Utils::sci(1.0 / psi_norm) << endl;
  }
  asymVecPtr->scale(1.0 / psi_norm);

  if (perturbSoln) {
    if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "\tIn LOCA::Bifurcation::PitchforkBord::ExtendedGroup::init(), applying random perturbation to initial solution of size:" 
	 << LOCA::Utils::sci(perturbSize) << endl;
    }
    NOX::Abstract::Vector *perturb = pfXVec.getXVec().clone(NOX::ShapeCopy);
    perturb->random();
    perturb->scale(pfXVec.getXVec());
    pfXVec.getXVec().update(perturbSize, *perturb, 1.0);
    grpPtr->setX(pfXVec.getXVec());
    delete perturb;
  }
}

double
LOCA::Bifurcation::PitchforkBord::ExtendedGroup::lTransNorm(
					const NOX::Abstract::Vector& n) const
{
  return lengthVecPtr->innerProduct(n) / lengthVecPtr->length();
}
