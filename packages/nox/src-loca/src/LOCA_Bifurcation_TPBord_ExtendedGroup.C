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

#include "LOCA_Bifurcation_TPBord_ExtendedGroup.H"
#include "LOCA_Bifurcation_TPBord_AbstractGroup.H"
#include "LOCA_Parameter_Vector.H"
#include "NOX_Parameter_List.H"
#include "LOCA_Utils.H"
#include "LOCA_ErrorCheck.H"

LOCA::Bifurcation::TPBord::ExtendedGroup::ExtendedGroup(
			      LOCA::Bifurcation::TPBord::AbstractGroup& g,
			      NOX::Parameter::List& bifParamList)
  : grpPtr(&g),
    tpXVec(g.getX(), g.getX(), 0.0),
    tpFVec(g.getX(), g.getX(), 0.0),
    tpNewtonVec(g.getX(), g.getX(), 0.0),
    lengthVecPtr(NULL), 
    bifParamId(0), 
    derivResidualParamPtr(NULL), 
    derivNullResidualParamPtr(NULL), 
    ownsGroup(false),
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false)
{
  const char *func = "LOCA::Bifurcation::TPBord::ExtendedGroup()";

  if (!bifParamList.isParameter("Bifurcation Parameter")) {
    LOCA::ErrorCheck::throwError(func,
				 "\"Bifurcation Parameter\" name is not set!");
  }
  string bifParamName = bifParamList.getParameter("Bifurcation Parameter",
						  "None");
  const ParameterVector& p = grpPtr->getParams();
  bifParamId = p.getIndex(bifParamName);

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

  lengthVecPtr = lenVecPtr->clone(NOX::DeepCopy);
  derivResidualParamPtr = lenVecPtr->clone(NOX::ShapeCopy);
  derivNullResidualParamPtr = lenVecPtr->clone(NOX::ShapeCopy);
  tpXVec.getNullVec() = *nullVecPtr;

  init(perturbSoln, perturbSize);
}

LOCA::Bifurcation::TPBord::ExtendedGroup::ExtendedGroup(
			      LOCA::Bifurcation::TPBord::AbstractGroup& g,
			      const NOX::Abstract::Vector& lenVec,
			      const NOX::Abstract::Vector& nullVec,
			      int paramId)
  : grpPtr(&g),
    tpXVec(g.getX(), nullVec, 0.0),
    tpFVec(lenVec, lenVec, 0.0),
    tpNewtonVec(lenVec, lenVec, 0.0),
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

LOCA::Bifurcation::TPBord::ExtendedGroup::ExtendedGroup(
			  const LOCA::Bifurcation::TPBord::AbstractGroup& g,
			  const NOX::Abstract::Vector& lenVec,
			  const NOX::Abstract::Vector& nullVec,
			  int paramId)
  : grpPtr(dynamic_cast<LOCA::Bifurcation::TPBord::AbstractGroup*>(g.clone())),
    tpXVec(g.getX(), lenVec, 0.0),
    tpFVec(lenVec, lenVec, 0.0),
    tpNewtonVec(lenVec, lenVec, 0.0),
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

LOCA::Bifurcation::TPBord::ExtendedGroup::ExtendedGroup(
		      const LOCA::Bifurcation::TPBord::ExtendedGroup& source, 
		      NOX::CopyType type)
  : grpPtr(dynamic_cast<LOCA::Bifurcation::TPBord::AbstractGroup*>(source.grpPtr->clone())), 
    tpXVec(source.tpXVec, type),
    tpFVec(source.tpFVec, type),
    tpNewtonVec(source.tpNewtonVec, type),
    lengthVecPtr(source.lengthVecPtr->clone(type)), 
    bifParamId(source.bifParamId),
    derivResidualParamPtr(source.derivResidualParamPtr->clone(type)),
    derivNullResidualParamPtr(source.derivNullResidualParamPtr->clone(type)),
    ownsGroup(true),
    isValidF(source.isValidF),
    isValidJacobian(source.isValidJacobian),
    isValidNewton(source.isValidNewton) {}


LOCA::Bifurcation::TPBord::ExtendedGroup::~ExtendedGroup() 
{
  if (ownsGroup)
    delete grpPtr;
  delete lengthVecPtr;
  delete derivResidualParamPtr;
  delete derivNullResidualParamPtr;
}

NOX::Abstract::Group&
LOCA::Bifurcation::TPBord::ExtendedGroup::operator=(const NOX::Abstract::Group& source)
{
  return *this = 
    dynamic_cast<const LOCA::Bifurcation::TPBord::ExtendedGroup&>(source);
}

LOCA::Continuation::AbstractGroup&
LOCA::Bifurcation::TPBord::ExtendedGroup::operator=(
			const LOCA::Continuation::AbstractGroup& source)
{
  return *this = 
    dynamic_cast<const LOCA::Bifurcation::TPBord::ExtendedGroup&>(source);
}

LOCA::Extended::AbstractGroup&
LOCA::Bifurcation::TPBord::ExtendedGroup::operator=(
			const LOCA::Extended::AbstractGroup& source)
{
  return *this = 
    dynamic_cast<const LOCA::Bifurcation::TPBord::ExtendedGroup&>(source);
}

LOCA::Bifurcation::TPBord::ExtendedGroup&
LOCA::Bifurcation::TPBord::ExtendedGroup::operator=(
		      const LOCA::Bifurcation::TPBord::ExtendedGroup& source) 
{

  // Protect against A = A
  if (this != &source) {
    
    // Copy values
    *grpPtr = *source.grpPtr;
    tpXVec = source.tpXVec;
    tpFVec = source.tpFVec;
    tpNewtonVec = source.tpNewtonVec;
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
LOCA::Bifurcation::TPBord::ExtendedGroup::clone(NOX::CopyType type) const 
{
  return new LOCA::Bifurcation::TPBord::ExtendedGroup(*this, type);
}

void
LOCA::Bifurcation::TPBord::ExtendedGroup::setParams(
					      const LOCA::ParameterVector& p) 
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  grpPtr->setParams(p);
}

const LOCA::ParameterVector&
LOCA::Bifurcation::TPBord::ExtendedGroup::getParams() const 
{
  return grpPtr->getParams();
}

void
LOCA::Bifurcation::TPBord::ExtendedGroup::setParam(int paramID, double val)
{
  grpPtr->setParam(paramID, val);
}

double
LOCA::Bifurcation::TPBord::ExtendedGroup::getParam(int paramID) const
{
  return grpPtr->getParam(paramID);
}

void
LOCA::Bifurcation::TPBord::ExtendedGroup::setParam(string paramID, double val)
{
  grpPtr->setParam(paramID, val);
}

double
LOCA::Bifurcation::TPBord::ExtendedGroup::getParam(string paramID) const
{
  return grpPtr->getParam(paramID);
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBord::ExtendedGroup::computeDfDp(int paramID, 
					      NOX::Abstract::Vector& result)
{
  string callingFunction = 
    "LOCA::Bifurcation::TPBord::ExtendedGroup::computeDfDp()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Cast result to TP vector
  LOCA::Bifurcation::TPBord::ExtendedVector& tp_result = 
    dynamic_cast<LOCA::Bifurcation::TPBord::ExtendedVector&>(result);

  // Make sure F is valid
  if (!isF()) {
    status = computeF();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Compute df/dp
  status = grpPtr->computeDfDp(paramID, tp_result.getXVec());
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute d(Jn)/dp
  status = grpPtr->computeDJnDp(tpXVec.getNullVec(), paramID,
				tpFVec.getNullVec(), tp_result.getNullVec());
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Set parameter componenet
  tp_result.getBifParam() = 0.0;

  return finalStatus;
}

double
LOCA::Bifurcation::TPBord::ExtendedGroup::getBifParam() const 
{
  return grpPtr->getParam(bifParamId);
}

void
LOCA::Bifurcation::TPBord::ExtendedGroup::setBifParam(double param) 
{
  grpPtr->setParam(bifParamId, param);

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

void
LOCA::Bifurcation::TPBord::ExtendedGroup::setX(const NOX::Abstract::Vector& y) 
{
  setX( dynamic_cast<const LOCA::Bifurcation::TPBord::ExtendedVector&>(y) );
}

void
LOCA::Bifurcation::TPBord::ExtendedGroup::setX(
			 const LOCA::Bifurcation::TPBord::ExtendedVector& y) 
{
  grpPtr->setX( y.getXVec() );
  tpXVec = y;
  setBifParam(tpXVec.getBifParam());

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

void
LOCA::Bifurcation::TPBord::ExtendedGroup::computeX(
					      const NOX::Abstract::Group& g, 
					      const NOX::Abstract::Vector& d,
					      double step) 
{
  computeX( dynamic_cast<const LOCA::Bifurcation::TPBord::ExtendedGroup&>(g),
	    dynamic_cast<const LOCA::Bifurcation::TPBord::ExtendedVector&>(d),
	    step );
}

void
LOCA::Bifurcation::TPBord::ExtendedGroup::computeX(
			    const LOCA::Bifurcation::TPBord::ExtendedGroup& g, 
			    const LOCA::Bifurcation::TPBord::ExtendedVector& d,
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
LOCA::Bifurcation::TPBord::ExtendedGroup::computeF() 
{
  if (isValidF)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::Bifurcation::TPBord::ExtendedGroup::computeF()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Compute underlying F
  if (!grpPtr->isF()) {
    status = grpPtr->computeF();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }
  tpFVec.getXVec() = grpPtr->getF();
  
  // Compute underlying Jacobian
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Compute J*n
  status = grpPtr->applyJacobian(tpXVec.getNullVec(), 
				 tpFVec.getNullVec());
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  
  // Compute phi^T*n
  tpFVec.getBifParam() = lTransNorm(tpXVec.getNullVec()) - 1.0;
  
  isValidF = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBord::ExtendedGroup::computeJacobian() 
{
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::Bifurcation::TPBord::ExtendedGroup::computeJacobian()";
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
  status = grpPtr->computeDJnDp(tpXVec.getNullVec(), 
				bifParamId,
				tpFVec.getNullVec(), 
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
LOCA::Bifurcation::TPBord::ExtendedGroup::computeGradient() 
{
  return NOX::Abstract::Group::NotDefined;
}
   
NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBord::ExtendedGroup::computeNewton(
						 NOX::Parameter::List& params) 
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::Bifurcation::TPBord::ExtendedGroup::computeNewton()";
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
  tpNewtonVec.init(0.0);

  status = applyJacobianInverse(params, tpFVec, tpNewtonVec);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  tpNewtonVec.scale(-1.0);
  isValidNewton = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBord::ExtendedGroup::applyJacobian(
				     const NOX::Abstract::Vector& input,
				     NOX::Abstract::Vector& result) const 
{
  string callingFunction = 
    "LOCA::Bifurcation::TPBord::ExtendedGroup::applyJacobian()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  if (!isJacobian()) {
    LOCA::ErrorCheck::throwError(callingFunction,
				 "Called with invalid Jacobian!");
  }

  // Cast vectors to TPBordVectors
  const LOCA::Bifurcation::TPBord::ExtendedVector& tp_input = 
    dynamic_cast<const LOCA::Bifurcation::TPBord::ExtendedVector&>(input);
  LOCA::Bifurcation::TPBord::ExtendedVector& tp_result = 
    dynamic_cast<LOCA::Bifurcation::TPBord::ExtendedVector&>(result);

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

  // compute J*x + p*dR/dp
  result_x.update(input_param, *derivResidualParamPtr, 1.0);

  // compute J*y
  status = grpPtr->applyJacobian(input_null, result_null);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // compute J*y + p*dJy/dp
  result_null.update(input_param, *derivNullResidualParamPtr, 1.0);

  // compute (dJy/dx)*x
  status = grpPtr->computeDJnDxa(tpXVec.getNullVec(), input_x, 
				 tpFVec.getNullVec(), *tmp);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  
  // compute (dJy/dx)*x + J*y + p*dJy/dp
  result_null.update(1.0, *tmp, 1.0);

  // compute l^T*y
  result_param = lTransNorm(input_null);

  delete tmp;

  return finalStatus;
}

// NOX::Abstract::Group::ReturnType
// LOCA::Bifurcation::TPBord::ExtendedGroup::applyJacobianMultiVector(
// 				     const NOX::Abstract::MultiVector& input,
// 				     NOX::Abstract::MultiVector& result) const 
// {
//   string callingFunction = 
//     "LOCA::Bifurcation::TPBord::ExtendedGroup::applyJacobianMultiVector()";
//   NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
//   NOX::Abstract::Group::ReturnType status;

//   if (!isJacobian()) {
//     LOCA::ErrorCheck::throwError(callingFunction,
// 				 "Called with invalid Jacobian!");
//   }

//   // Cast vectors to TPBordVectors
//   const LOCA::Bifurcation::TPBord::ExtendedMultiVector& tp_input = 
//     dynamic_cast<const LOCA::Bifurcation::TPBord::ExtendedMultiVector&>(input);
//   LOCA::Bifurcation::TPBord::ExtendedMultiVector& tp_result = 
//     dynamic_cast<LOCA::Bifurcation::TPBord::ExtendedMultiVector&>(result);

//   // Get constant references to input vector components
//   const NOX::Abstract::MultiVector& input_x = tp_input.getXMultiVec();
//   const NOX::Abstract::MultiVector& input_null = tp_input.getNullMultiVec();
//   const NOX::Abstract::MultiVector::DenseMatrix& input_param = 
//     tp_input.getBifParamMatrix();

//   // Get non-constant references to result vector components
//   NOX::Abstract::MultiVector& result_x = tp_result.getXMultiVec();
//   NOX::Abstract::MultiVector& result_null = tp_result.getNullMultiVec();
//   NOX::Abstract::MultiVector::DenseMatrix& result_param = 
//     tp_result.getBifParamMatrix();

//   // Temporary vector
//   NOX::Abstract::MultiVector *tmp = input_null.clone(NOX::ShapeCopy);

//   // verify underlying Jacobian is valid
//   if (!grpPtr->isJacobian()) {
//     status = grpPtr->computeJacobian();
//     finalStatus = 
//       LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
// 						   callingFunction);
//   }

//   // compute J*x
//   status = grpPtr->applyJacobianMultiVector(input_x, result_x);
//   finalStatus = 
//     LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
// 						 callingFunction);

//   // compute J*x + p*dR/dp
//   result_x.update(1.0, *derivResidualParamPtr, input_param, 1.0);

//   // compute J*y
//   status = grpPtr->applyJacobianMulti(input_null, result_null);
//   finalStatus = 
//     LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
// 						 callingFunction);

//   // compute J*y + p*dJy/dp
//   result_null.update(1.0, *derivNullResidualParamPtr, input_param, 1.0);

//   // compute (dJy/dx)*x
//   status = grpPtr->computeDJnDxa(tpXVec.getNullVec(), input_x, 
// 				 tpFVec.getNullVec(), *tmp);
//   finalStatus = 
//     LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
// 						 callingFunction);
  
//   // compute (dJy/dx)*x + J*y + p*dJy/dp
//   result_null.update(1.0, *tmp, 1.0);

//   // compute l^T*y
//   result_param = lTransNorm(input_null);

//   delete tmp;

//   return finalStatus;
// }

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBord::ExtendedGroup::applyJacobianTranspose(
					  const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBord::ExtendedGroup::applyJacobianInverse(
					NOX::Parameter::List& params,
					const NOX::Abstract::Vector& input,
					NOX::Abstract::Vector& result) const 
{
  const NOX::Abstract::Vector* inputs[1];
  NOX::Abstract::Vector* results[1];

  inputs[0] = &input;
  results[0] = &result;

  return  applyJacobianInverseMulti(params, inputs, results, 1);
}

NOX::Abstract::Group::ReturnType 
LOCA::Bifurcation::TPBord::ExtendedGroup::applyRightPreconditioning(
					 bool useTranspose,
					 NOX::Parameter::List& params,
					 const NOX::Abstract::Vector& input, 
					 NOX::Abstract::Vector& result) const
{
  if (useTranspose) {
    LOCA::ErrorCheck::printWarning(
	"LOCA::Bifurcation::TPBord::ExtendedGroup::applyRightPreconditioning",
	"Transpose of right preconditioner not implemented");
    return NOX::Abstract::Group::NotDefined;
  }

  string callingFunction = 
    "LOCA::Bifurcation::TPBord::ExtendedGroup::applyRightPreconditioning()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  if (!isJacobian()) {
    LOCA::ErrorCheck::throwError(callingFunction,
				 "Called with invalid Jacobian!");
  }

  // cast vectors to turning point vectors
  const LOCA::Bifurcation::TPBord::ExtendedVector& tp_input = 
    dynamic_cast<const LOCA::Bifurcation::TPBord::ExtendedVector&>(input);
  LOCA::Bifurcation::TPBord::ExtendedVector& tp_result = 
    dynamic_cast<LOCA::Bifurcation::TPBord::ExtendedVector&>(result);

  // Get componenets of input vector
  const NOX::Abstract::Vector& input_x = tp_input.getXVec();
  const NOX::Abstract::Vector& input_y = tp_input.getNullVec();
  double input_p = tp_input.getBifParam();

  // Get components of result vector.  Note these are references.
  NOX::Abstract::Vector& result_x = tp_result.getXVec();
  NOX::Abstract::Vector& result_y = tp_result.getNullVec();
  double& result_p = tp_result.getBifParam();

  // Temporary vectors
  NOX::Abstract::Vector* a = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector* b = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector* c = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector* d = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector* tmp = input_x.clone(NOX::ShapeCopy);

  // Solve P*a = input_x
  finalStatus = grpPtr->applyRightPreconditioning(useTranspose, params, 
						  input_x, *a);
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  // Solve P*b = dF/dp
  status = grpPtr->applyRightPreconditioning(useTranspose, params, 
					     *derivResidualParamPtr, *b);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  
  // Compute (dJy/dx)*a - input_y
  status = grpPtr->computeDJnDxa(tpXVec.getNullVec(), *a, tpFVec.getNullVec(), 
				 *tmp);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  tmp->update(-1.0, input_y, 1.0);

  // Solve P*c = (dJy/dx)*a - input_y
  status = grpPtr->applyRightPreconditioning(useTranspose, params, *tmp, *c);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute (dJy/dx)*b - dJy/dp
  status = grpPtr->computeDJnDxa(tpXVec.getNullVec(), *b, tpFVec.getNullVec(),
				 *tmp);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  tmp->update(-1.0, *derivNullResidualParamPtr, 1.0);

  // Solve P*d = (dJy/dx)*b - dJy/dp
  status = grpPtr->applyRightPreconditioning(useTranspose, params, *tmp, *d);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  double w = (input_p + lTransNorm(*c))/lTransNorm(*d);
  
  result_x.update(1.0, *a, -w, *b, 0.0);
  result_y.update(-1.0, *c, w, *d, 0.0);
  result_p = w;

  delete a;
  delete b;
  delete c;
  delete d;
  delete tmp;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBord::ExtendedGroup::applyJacobianInverseMulti(
			    NOX::Parameter::List& params,
			    const NOX::Abstract::Vector* const* inputs,
			    NOX::Abstract::Vector** results, int nVecs) const 
{

  string callingFunction = 
    "LOCA::Bifurcation::TPBord::ExtendedGroup::applyJacobianInverseMulti()";
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
    new const NOX::Abstract::Vector*[m+1];
  const NOX::Abstract::Vector** inputs_null =
    new const NOX::Abstract::Vector*[m+1];
  double* inputs_params = new double[m];

  NOX::Abstract::Vector** tmp1 = new NOX::Abstract::Vector*[m+1];
  NOX::Abstract::Vector** tmp2 = new NOX::Abstract::Vector*[m+1];
  NOX::Abstract::Vector** tmp3 = new NOX::Abstract::Vector*[m+1];

  const LOCA::Bifurcation::TPBord::ExtendedVector* constTPVecPtr;

  for (int i=0; i<m; i++) {
    constTPVecPtr = 
      dynamic_cast<const LOCA::Bifurcation::TPBord::ExtendedVector*>(inputs[i]);
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

  // verify underlying Jacobian is valid
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Solve J*tmp1 = inputs_x
  status = grpPtr->applySingularJacobianInverseMulti(params, inputs_x, 
						     tpXVec.getNullVec(),
						     tpFVec.getNullVec(),
						     tmp1, m+1);
  finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);

  // Compute tmp2 = (dJy/dx)*tmp1 - inputs_null
  for (int i=0; i<m+1; i++) {
    status = grpPtr->computeDJnDxa(tpXVec.getNullVec(), *tmp1[i],
				   tpFVec.getNullVec(), *tmp2[i]);

    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);

    tmp2[i]->update(-1.0, *inputs_null[i], 1.0);
  } 

  // verify underlying Jacobian is valid
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Solve J*tmp3 = tmp2
  status = grpPtr->applySingularJacobianInverseMulti(params, tmp2, 
						     tpXVec.getNullVec(),
						     tpFVec.getNullVec(),
						     tmp3, m+1);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  
  // Compute and set results
  double denom = lTransNorm(*tmp3[m]);
  double w;
  LOCA::Bifurcation::TPBord::ExtendedVector* tpVecPtr;
  for (int i=0; i<m; i++) {
    tpVecPtr = 
      dynamic_cast<LOCA::Bifurcation::TPBord::ExtendedVector*>(results[i]);

    w = (inputs_params[i] + lTransNorm(*tmp3[i]))/denom;
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

  return finalStatus;
}

bool
LOCA::Bifurcation::TPBord::ExtendedGroup::isF() const 
{
  return isValidF;
}

bool
LOCA::Bifurcation::TPBord::ExtendedGroup::isJacobian() const 
{
  return isValidJacobian;
}

bool
LOCA::Bifurcation::TPBord::ExtendedGroup::isGradient() const 
{
  return false;
}

bool
LOCA::Bifurcation::TPBord::ExtendedGroup::isNewton() const 
{
  return isValidNewton;
}
  
const NOX::Abstract::Vector&
LOCA::Bifurcation::TPBord::ExtendedGroup::getX() const 
{
  return tpXVec;
}

const NOX::Abstract::Vector&
LOCA::Bifurcation::TPBord::ExtendedGroup::getF() const 
{
  return tpFVec;
}

double
LOCA::Bifurcation::TPBord::ExtendedGroup::getNormF() const 
{
  return tpFVec.norm();
}

const NOX::Abstract::Vector&
LOCA::Bifurcation::TPBord::ExtendedGroup::getGradient() const 
{
  LOCA::ErrorCheck::throwError(
		   "LOCA::Bifurcation::TPBord::ExtendedGroup::getGradient()",
		   " - not implemented");
  return getNewton();
}

const NOX::Abstract::Vector&
LOCA::Bifurcation::TPBord::ExtendedGroup::getNewton() const 
{
  return tpNewtonVec;
}

double
LOCA::Bifurcation::TPBord::ExtendedGroup::getNormNewtonSolveResidual() const 
{
  string callingFunction = 
    "LOCA::Bifurcation::TPBord::ExtendedGroup::getNormNewtonSolveResidual()";
  NOX::Abstract::Group::ReturnType finalStatus;
  LOCA::Bifurcation::TPBord::ExtendedVector residual = tpFVec;
  
  finalStatus = applyJacobian(tpNewtonVec, residual);
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  residual.update(1.0, tpFVec, 1.0);
  return residual.norm();
}

void
LOCA::Bifurcation::TPBord::ExtendedGroup::printSolution(const double conParam) const 
{
  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "LOCA::Bifurcation::TPBord::ExtendedGroup::printSolution\n";

    cout << "Turning Point located at: " << LOCA::Utils::sci(conParam) 
	 << "   " 
	 << LOCA::Utils::sci(getBifParam())
	 << endl;

    cout << "\tPrinting Solution Vector for conParam = " 
	 << LOCA::Utils::sci(conParam) << endl;
  }
  grpPtr->printSolution(conParam);
  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "\tPrinting Null Vector for bif param = " 
	 << LOCA::Utils::sci(getBifParam()) << endl;
  }
  grpPtr->printSolution(tpXVec.getNullVec(), tpXVec.getBifParam());
}

const LOCA::Continuation::AbstractGroup&
LOCA::Bifurcation::TPBord::ExtendedGroup::getUnderlyingGroup() const
{
  return *grpPtr;
}

LOCA::Continuation::AbstractGroup&
LOCA::Bifurcation::TPBord::ExtendedGroup::getUnderlyingGroup()
{
  return *grpPtr;
}

void
LOCA::Bifurcation::TPBord::ExtendedGroup::init(bool perturbSoln, 
					       double perturbSize)
{
  tpXVec.getBifParam() = getBifParam();

  // Rescale length vector so that the normalization condition is met
  double lVecDotNullVec = lTransNorm(tpXVec.getNullVec());

  if (lVecDotNullVec == 0.0) {
    LOCA::ErrorCheck::throwError(
		   "LOCA::Bifurcation::TPBord::ExtendedGroup::init()",
		   "null vector can be orthogonal to length-scaling vector");
  }
  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "\tIn LOCA::Bifurcation::TPBord::ExtendedGroup::init(), scaling null vector by:" 
	 << LOCA::Utils::sci(1.0 / lVecDotNullVec) << endl;
  }
  tpXVec.getNullVec().scale(1.0/lVecDotNullVec);

  if (perturbSoln) {
    if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "\tIn LOCA::Bifurcation::TPBord::ExtendedGroup::init(), applying random perturbation to initial solution of size:" 
	 << LOCA::Utils::sci(perturbSize) << endl;
    }
    NOX::Abstract::Vector *perturb = tpXVec.getXVec().clone(NOX::ShapeCopy);
    perturb->random();
    perturb->scale(tpXVec.getXVec());
    tpXVec.getXVec().update(perturbSize, *perturb, 1.0);
    delete perturb;
  }
}

double
LOCA::Bifurcation::TPBord::ExtendedGroup::lTransNorm(
					const NOX::Abstract::Vector& n) const
{
  return lengthVecPtr->dot(n) / lengthVecPtr->length();
}
