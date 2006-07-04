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

#include "LOCA_Bifurcation_HopfBord_ExtendedGroup.H"
#include "LOCA_Bifurcation_HopfBord_AbstractGroup.H" 
#include "LOCA_Parameter_Vector.H"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Utils.H"

LOCA::Bifurcation::HopfBord::ExtendedGroup::ExtendedGroup(
			      LOCA::Bifurcation::HopfBord::AbstractGroup& g,
			      Teuchos::ParameterList& bifParamList)
  : grpPtr(&g),
    hopfXVec(g.getX(), g.getX(), g.getX(), 0.0, 0.0),
    hopfFVec(g.getX(), g.getX(), g.getX(), 0.0, 0.0),
    hopfNewtonVec(g.getX(), g.getX(), g.getX(), 0.0, 0.0),
    lengthVecPtr(NULL), 
    bifParamId(0), 
    derivResidualParamPtr(NULL), 
    derivRealEigenResidualParamPtr(NULL),
    derivImagEigenResidualParamPtr(NULL),
    massTimesYPtr(NULL),
    minusMassTimesZPtr(NULL),
    ownsGroup(false),
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false)
{
  const char *func = "LOCA::Bifurcation::HopfBord::ExtendedGroup()";

  if (!bifParamList.isParameter("Bifurcation Parameter")) {
    LOCA::ErrorCheck::throwError(func,
				 "\"Bifurcation Parameter\" name is not set!");
  }
  string bifParamName = bifParamList.get("Bifurcation Parameter",
						  "None");
  const ParameterVector& p = grpPtr->getParams();
  bifParamId = p.getIndex(bifParamName);

  if (!bifParamList.isParameter("Length Normalization Vector")) {
    LOCA::ErrorCheck::throwError(func,
			   "\"Length Normalization Vector\" is not set!");
  }
  NOX::Abstract::Vector* lenVecPtr = 
    bifParamList.getAnyPtrParameter<NOX::Abstract::Vector>("Length Normalization Vector");

  if (!bifParamList.isParameter("Initial Real Eigenvector")) {
    LOCA::ErrorCheck::throwError(func,
			   "\"Initial Real Eigenvector\" is not set!");
  }
  const NOX::Abstract::Vector* realEigenVecPtr = 
    bifParamList.getAnyConstPtrParameter<NOX::Abstract::Vector>("Initial Real Eigenvector");

  if (!bifParamList.isParameter("Initial Imaginary Eigenvector")) {
    LOCA::ErrorCheck::throwError(func,
			   "\"Initial Imaginary Eigenvector\" is not set!");
  }
  const NOX::Abstract::Vector* imagEigenVecPtr = 
    bifParamList.getAnyConstPtrParameter<NOX::Abstract::Vector>("Initial Imaginary Eigenvector");

  if (!bifParamList.isParameter("Initial Frequency")) {
    LOCA::ErrorCheck::throwError(func,
				 "\"Initial Frequency\" name is not set!");
  }
  double frequency = bifParamList.get("Initial Frequency",0.0);

  bool perturbSoln = bifParamList.get("Perturb Initial Solution", 
					       false);
  double perturbSize = bifParamList.get("Relative Perturbation Size", 
						 1.0e-3);

  hopfXVec.getRealEigenVec() = *realEigenVecPtr;
  hopfXVec.getImagEigenVec() = *imagEigenVecPtr;
  hopfXVec.getFrequency() = frequency;
  lengthVecPtr = lenVecPtr;
  derivResidualParamPtr = lenVecPtr->clone(NOX::ShapeCopy);
  derivRealEigenResidualParamPtr = lenVecPtr->clone(NOX::ShapeCopy);
  derivImagEigenResidualParamPtr = lenVecPtr->clone(NOX::ShapeCopy);
  massTimesYPtr = lenVecPtr->clone(NOX::ShapeCopy);
  minusMassTimesZPtr = lenVecPtr->clone(NOX::ShapeCopy);

  init(perturbSoln, perturbSize);
}

LOCA::Bifurcation::HopfBord::ExtendedGroup::ExtendedGroup(
			      LOCA::Bifurcation::HopfBord::AbstractGroup& g,
			      const NOX::Abstract::Vector& realEigenVec,
			      const NOX::Abstract::Vector& imaginaryEigenVec,
			      NOX::Abstract::Vector& lenVec,
			      double frequency,
			      int paramId)
  : grpPtr(&g), 
    hopfXVec(g.getX(), realEigenVec, imaginaryEigenVec, frequency, 0.0),
    hopfFVec(g.getX(), realEigenVec, imaginaryEigenVec, frequency, 0.0),
    hopfNewtonVec(g.getX(), realEigenVec, imaginaryEigenVec, frequency, 0.0),
    lengthVecPtr(&lenVec), 
    bifParamId(paramId), 
    derivResidualParamPtr(lenVec.clone(NOX::ShapeCopy)), 
    derivRealEigenResidualParamPtr(lenVec.clone(NOX::ShapeCopy)),
    derivImagEigenResidualParamPtr(lenVec.clone(NOX::ShapeCopy)),
    massTimesYPtr(lenVec.clone(NOX::ShapeCopy)),
    minusMassTimesZPtr(lenVec.clone(NOX::ShapeCopy)),
    ownsGroup(false),
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false)
{
  init();
}

LOCA::Bifurcation::HopfBord::ExtendedGroup::ExtendedGroup(
			   const LOCA::Bifurcation::HopfBord::AbstractGroup& g,
			   const NOX::Abstract::Vector& realEigenVec,
			   const NOX::Abstract::Vector& imaginaryEigenVec,
			   NOX::Abstract::Vector& lenVec,
			   double frequency,
			   int paramId)
  : grpPtr(dynamic_cast<LOCA::Bifurcation::HopfBord::AbstractGroup*>(g.clone())), 
    hopfXVec(g.getX(), realEigenVec, imaginaryEigenVec, frequency, 0.0),
    hopfFVec(g.getX(), realEigenVec, imaginaryEigenVec, frequency, 0.0),
    hopfNewtonVec(g.getX(), realEigenVec, imaginaryEigenVec, frequency, 0.0),
    lengthVecPtr(&lenVec), 
    bifParamId(paramId), 
    derivResidualParamPtr(lenVec.clone(NOX::ShapeCopy)), 
    derivRealEigenResidualParamPtr(lenVec.clone(NOX::ShapeCopy)),
    derivImagEigenResidualParamPtr(lenVec.clone(NOX::ShapeCopy)),
    massTimesYPtr(lenVec.clone(NOX::ShapeCopy)),
    minusMassTimesZPtr(lenVec.clone(NOX::ShapeCopy)),
    ownsGroup(true),
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false)
{
  init();
}

LOCA::Bifurcation::HopfBord::ExtendedGroup::ExtendedGroup(
		    const LOCA::Bifurcation::HopfBord::ExtendedGroup& source, 
		    NOX::CopyType type)
  : grpPtr(dynamic_cast<LOCA::Bifurcation::HopfBord::AbstractGroup*>(source.grpPtr->clone())), 
    hopfXVec(source.hopfXVec, type),
    hopfFVec(source.hopfFVec, type),
    hopfNewtonVec(source.hopfNewtonVec, type),
    lengthVecPtr(source.lengthVecPtr), 
    bifParamId(source.bifParamId),
    derivResidualParamPtr(source.derivResidualParamPtr->clone(type)),
    derivRealEigenResidualParamPtr(source.derivRealEigenResidualParamPtr->clone(type)),
    derivImagEigenResidualParamPtr(source.derivImagEigenResidualParamPtr->clone(type)),
    massTimesYPtr(source.massTimesYPtr->clone(type)),
    minusMassTimesZPtr(source.minusMassTimesZPtr->clone(type)),
    ownsGroup(true),
    isValidF(source.isValidF),
    isValidJacobian(source.isValidJacobian),
    isValidNewton(source.isValidNewton) {}


LOCA::Bifurcation::HopfBord::ExtendedGroup::~ExtendedGroup() 
{
  if (ownsGroup)
    delete grpPtr;
  delete derivResidualParamPtr;
  delete derivRealEigenResidualParamPtr;
  delete derivImagEigenResidualParamPtr;
  delete massTimesYPtr;
  delete minusMassTimesZPtr;
}

LOCA::Continuation::AbstractGroup&
LOCA::Bifurcation::HopfBord::ExtendedGroup::operator=(
			     const LOCA::Continuation::AbstractGroup& source)
{
  return *this = 
    dynamic_cast<const LOCA::Bifurcation::HopfBord::ExtendedGroup&>(source);
}

LOCA::Extended::AbstractGroup&
LOCA::Bifurcation::HopfBord::ExtendedGroup::operator=(
			     const LOCA::Extended::AbstractGroup& source)
{
  return *this = 
    dynamic_cast<const LOCA::Bifurcation::HopfBord::ExtendedGroup&>(source);
}

NOX::Abstract::Group&
LOCA::Bifurcation::HopfBord::ExtendedGroup::operator=(
				         const NOX::Abstract::Group& source)
{
  return *this = 
    dynamic_cast<const LOCA::Bifurcation::HopfBord::ExtendedGroup&>(source);
}

LOCA::Bifurcation::HopfBord::ExtendedGroup&
LOCA::Bifurcation::HopfBord::ExtendedGroup::operator=(
		    const LOCA::Bifurcation::HopfBord::ExtendedGroup& source) 
{

  // Protect against A = A
  if (this != &source) {

    *grpPtr = *source.grpPtr;
    hopfXVec = source.hopfXVec;
    hopfFVec = source.hopfFVec;
    hopfNewtonVec = source.hopfNewtonVec;
    *lengthVecPtr = *source.lengthVecPtr;
    bifParamId = source.bifParamId;
    *derivResidualParamPtr = *source.derivResidualParamPtr;
    *derivRealEigenResidualParamPtr = *source.derivRealEigenResidualParamPtr;
    *derivImagEigenResidualParamPtr = *source.derivImagEigenResidualParamPtr;
    *massTimesYPtr = *source.massTimesYPtr;
    *minusMassTimesZPtr = *source.minusMassTimesZPtr;
    ownsGroup = source.ownsGroup;
    isValidF = source.isValidF;
    isValidJacobian = source.isValidJacobian;
    isValidNewton = source.isValidNewton;
  }

  return *this;
}

NOX::Abstract::Group*
LOCA::Bifurcation::HopfBord::ExtendedGroup::clone(NOX::CopyType type) const 
{
  return new LOCA::Bifurcation::HopfBord::ExtendedGroup(*this, type);
}

void
LOCA::Bifurcation::HopfBord::ExtendedGroup::setParams(
					      const LOCA::ParameterVector& p) 
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  grpPtr->setParams(p);
}

const LOCA::ParameterVector&
LOCA::Bifurcation::HopfBord::ExtendedGroup::getParams() const 
{
  return grpPtr->getParams();
}

void
LOCA::Bifurcation::HopfBord::ExtendedGroup::setParam(int paramID, double val)
{
  grpPtr->setParam(paramID, val);
}

double
LOCA::Bifurcation::HopfBord::ExtendedGroup::getParam(int paramID) const
{
  return grpPtr->getParam(paramID);
}

void
LOCA::Bifurcation::HopfBord::ExtendedGroup::setParam(string paramID, 
						     double val)
{
  grpPtr->setParam(paramID, val);
}

double
LOCA::Bifurcation::HopfBord::ExtendedGroup::getParam(string paramID) const
{
  return grpPtr->getParam(paramID);
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::ExtendedGroup::computeDfDp(int paramID, 
					      NOX::Abstract::Vector& result)
{
  string callingFunction = 
    "LOCA::Bifurcation::HopfBord::ExtendedGroup::computeDfDp()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Cast result to Hopf vector
  LOCA::Bifurcation::HopfBord::ExtendedVector& hopf_result = 
    dynamic_cast<LOCA::Bifurcation::HopfBord::ExtendedVector&>(result);

  // Make sure F is valid
  if (!isF()) {
    status = computeF();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Compute df/dp
  status = grpPtr->computeDfDp(paramID, hopf_result.getXVec());
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute [(J+iwB)*(y+iz)]_p
  status = grpPtr->computeDCeDp(hopfXVec.getRealEigenVec(),
				hopfXVec.getImagEigenVec(),
				hopfXVec.getFrequency(),
				paramID,
				hopfFVec.getRealEigenVec(),
				hopfFVec.getImagEigenVec(),
				hopf_result.getRealEigenVec(),
				hopf_result.getImagEigenVec());
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Set frequency  component
  hopf_result.getFrequency() = 0.0;

  // Set parameter component
  hopf_result.getBifParam() = 0.0;

  return finalStatus;
}

double
LOCA::Bifurcation::HopfBord::ExtendedGroup::getBifParam() const 
{
  return grpPtr->getParam(bifParamId);
}

void
LOCA::Bifurcation::HopfBord::ExtendedGroup::setBifParam(double param) 
{
  grpPtr->setParam(bifParamId, param);

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

void
LOCA::Bifurcation::HopfBord::ExtendedGroup::setX(
					     const NOX::Abstract::Vector& y) 
{
  setX( dynamic_cast<const LOCA::Bifurcation::HopfBord::ExtendedVector&>(y) );
}

void
LOCA::Bifurcation::HopfBord::ExtendedGroup::setX(
			const LOCA::Bifurcation::HopfBord::ExtendedVector& y) 
{
  grpPtr->setX( y.getXVec() );
  hopfXVec = y;
  setBifParam(hopfXVec.getBifParam());

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

void
LOCA::Bifurcation::HopfBord::ExtendedGroup::computeX(
					     const NOX::Abstract::Group& g, 
					     const NOX::Abstract::Vector& d,
					     double step) 
{
  computeX(dynamic_cast<const LOCA::Bifurcation::HopfBord::ExtendedGroup&>(g),
	   dynamic_cast<const LOCA::Bifurcation::HopfBord::ExtendedVector&>(d),
	   step);
}

void
LOCA::Bifurcation::HopfBord::ExtendedGroup::computeX(
			 const LOCA::Bifurcation::HopfBord::ExtendedGroup& g, 
			 const LOCA::Bifurcation::HopfBord::ExtendedVector& d,
			 double step) 
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  grpPtr->computeX(*(g.grpPtr), d.getXVec(), step);
  hopfXVec.update(1.0, g.getX(), step, d, 0.0);
  setBifParam(hopfXVec.getBifParam());
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::ExtendedGroup::computeF() 
{
  if (isValidF)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::Bifurcation::HoopfBord::ExtendedGroup::computeF()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Compute underlying F
  if (!grpPtr->isF()) {
    status = grpPtr->computeF();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }
  hopfFVec.getXVec() = grpPtr->getF();
  
  // Compute underlying Jacobian
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Compute underlying Mass Matrix
  if (!grpPtr->isMassMatrix()) {
    status = grpPtr->computeMassMatrix();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Compute (J+iwB)*(y+iz)
  status = grpPtr->applyComplex(hopfXVec.getRealEigenVec(),
				hopfXVec.getImagEigenVec(),
				hopfXVec.getFrequency(),
				hopfFVec.getRealEigenVec(),
				hopfFVec.getImagEigenVec());
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  
  // Compute l^T*y - 1
  hopfFVec.getFrequency() = lTransNorm(hopfXVec.getRealEigenVec())-1.0;

  // Compute l^T*z
  hopfFVec.getBifParam() = lTransNorm(hopfXVec.getImagEigenVec());
  
  isValidF = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::ExtendedGroup::computeJacobian() 
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

  // Compute underlying [(J+iwB)*(y+iz)]_p
  status = grpPtr->computeDCeDp(hopfXVec.getRealEigenVec(),
				hopfXVec.getImagEigenVec(),
				hopfXVec.getFrequency(),
				bifParamId,
				hopfFVec.getRealEigenVec(),
				hopfFVec.getImagEigenVec(),
				*derivRealEigenResidualParamPtr,
				*derivImagEigenResidualParamPtr);
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

  // Compute underlying mass matrix
  if (!grpPtr->isMassMatrix()) {
    status = grpPtr->computeMassMatrix();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Compute B*y
  status = grpPtr->applyMassMatrix(hopfXVec.getRealEigenVec(),
				   *massTimesYPtr);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute -B*z
  status = grpPtr->applyMassMatrix(hopfXVec.getImagEigenVec(),
				   *minusMassTimesZPtr);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  minusMassTimesZPtr->scale(-1.0);

  isValidJacobian = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::ExtendedGroup::computeGradient() 
{
  return NOX::Abstract::Group::NotDefined;
}
   
NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::ExtendedGroup::computeNewton(
						Teuchos::ParameterList& params) 
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::Bifurcation::HopfBord::ExtendedGroup::computeNewton()";
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
  hopfNewtonVec.init(0.0);

  status = applyJacobianInverse(params, hopfFVec, hopfNewtonVec);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  hopfNewtonVec.scale(-1.0);
  isValidNewton = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::ExtendedGroup::applyJacobian(
					  const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  string callingFunction = 
    "LOCA::Bifurcation::HopfBord::ExtendedGroup::applyJacobian()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  if (!isJacobian()) {
    LOCA::ErrorCheck::throwError(callingFunction,
				 "Called with invalid Jacobian!");
  }

  // Cast vectors to HopfBordVectors
  const LOCA::Bifurcation::HopfBord::ExtendedVector& hopf_input = 
    dynamic_cast<const LOCA::Bifurcation::HopfBord::ExtendedVector&>(input);
  LOCA::Bifurcation::HopfBord::ExtendedVector& hopf_result = 
    dynamic_cast<LOCA::Bifurcation::HopfBord::ExtendedVector&>(result);

  // Get constant references to input vector components
  const NOX::Abstract::Vector& input_x = hopf_input.getXVec();
  const NOX::Abstract::Vector& input_y = hopf_input.getRealEigenVec();
  const NOX::Abstract::Vector& input_z = hopf_input.getImagEigenVec();
  double input_w = hopf_input.getFrequency();
  double input_param = hopf_input.getBifParam();

  // Get non-constant references to result vector components
  NOX::Abstract::Vector& result_x = hopf_result.getXVec();
  NOX::Abstract::Vector& result_y = hopf_result.getRealEigenVec();
  NOX::Abstract::Vector& result_z = hopf_result.getImagEigenVec();
  double& result_w = hopf_result.getFrequency();
  double& result_param = hopf_result.getBifParam();

  // Temporary vectors
  NOX::Abstract::Vector *tmp_real = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector *tmp_imag = input_x.clone(NOX::ShapeCopy);

  // verify underlying Jacobian is valid
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // verify underlying mass matrix is valid
  if (!grpPtr->isMassMatrix()) {
    status = grpPtr->computeMassMatrix();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // compute J*X
  finalStatus = grpPtr->applyJacobian(input_x, result_x);
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  // compute J*X + P*dR/dp
  result_x.update(input_param, *derivResidualParamPtr, 1.0);

  // compute (J+iwB)*(Y+iZ)
  status = grpPtr->applyComplex(input_y, input_z, hopfXVec.getFrequency(), 
				result_y, result_z);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // compute (J+iwB)*(Y+iZ) + P*d(J+iwB)*(y+iz)/dp + iW*B*(y+iz)
  result_y.update(input_param, *derivRealEigenResidualParamPtr, 
		  input_w, *minusMassTimesZPtr, 1.0);
  result_z.update(input_param, *derivImagEigenResidualParamPtr, 
		  input_w, *massTimesYPtr, 1.0);

  // compute (d(J+iwB)*(y+iz)/dx)*X
  status = grpPtr->computeDCeDxa(hopfXVec.getRealEigenVec(), 
				 hopfXVec.getImagEigenVec(), 
				 hopfXVec.getFrequency(),
				 input_x, 
				 hopfFVec.getRealEigenVec(), 
				 hopfFVec.getImagEigenVec(),
				 *tmp_real,
				 *tmp_imag);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  
  // (d(J+iwB)*(y+iz)/dx)*X + (J+iwB)*(Y+iZ) + P*d(J+iwB)*(y+iz)/dp + iW*B*(y+iz)
  result_y.update(1.0, *tmp_real, 1.0);
  result_z.update(1.0, *tmp_imag, 1.0);

  // compute phi^T*Y, phi^T*Z
  result_w = lTransNorm(input_y);
  result_param = lTransNorm(input_z);

  delete tmp_real;
  delete tmp_imag;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::ExtendedGroup::applyJacobianTranspose(
					 const NOX::Abstract::Vector& input,
					 NOX::Abstract::Vector& result) const 
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::ExtendedGroup::applyJacobianInverse(
					 Teuchos::ParameterList& params,
					 const NOX::Abstract::Vector& input,
					 NOX::Abstract::Vector& result) const 
{
  NOX::Abstract::Group::ReturnType res;
  const NOX::Abstract::Vector* inputs[1];
  NOX::Abstract::Vector* results[1];

  inputs[0] = &input;
  results[0] = &result;

  res = applyJacobianInverseMulti(params, inputs, results, 1);

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::ExtendedGroup::applyJacobianInverseMulti(
			    Teuchos::ParameterList& params,
			    const NOX::Abstract::Vector* const* inputs,
			    NOX::Abstract::Vector** results, int nVecs) const 
{
  string callingFunction = 
    "LOCA::Bifurcation::HopfBord::ExtendedGroup::applyJacobianInverseMulti()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  if (!isJacobian()) {
    LOCA::ErrorCheck::throwError(callingFunction,
				 "Called with invalid Jacobian!");
  }

  // Number of input vectors
  int m = nVecs; 
  
  // Build arrays of solution, eigenvector and parameter components
  const NOX::Abstract::Vector** inputs_x = 
    new const NOX::Abstract::Vector*[m+1];
  const NOX::Abstract::Vector** inputs_y =
    new const NOX::Abstract::Vector*[m+1];
  const NOX::Abstract::Vector** inputs_z =
    new const NOX::Abstract::Vector*[m+1];
  double* inputs_w = new double[m];
  double* inputs_params = new double[m];

  NOX::Abstract::Vector** results_x = new NOX::Abstract::Vector*[m+1];
  NOX::Abstract::Vector** results_y = new NOX::Abstract::Vector*[m+2];
  NOX::Abstract::Vector** results_z = new NOX::Abstract::Vector*[m+2];
  double** results_w = new double*[m];
  double** results_params = new double*[m];

  NOX::Abstract::Vector** tmp_real = new NOX::Abstract::Vector*[m+2];
  NOX::Abstract::Vector** tmp_imag = new NOX::Abstract::Vector*[m+2];

  const LOCA::Bifurcation::HopfBord::ExtendedVector* constHopfVecPtr;
  LOCA::Bifurcation::HopfBord::ExtendedVector* hopfVecPtr;

  for (int i=0; i<m; i++) {
    constHopfVecPtr = 
      dynamic_cast<const LOCA::Bifurcation::HopfBord::ExtendedVector*>(inputs[i]);
    inputs_x[i] = &(constHopfVecPtr->getXVec());
    inputs_y[i] = &(constHopfVecPtr->getRealEigenVec());
    inputs_z[i] = &(constHopfVecPtr->getImagEigenVec());
    inputs_w[i] = constHopfVecPtr->getFrequency();
    inputs_params[i] = constHopfVecPtr->getBifParam();
 
    hopfVecPtr = 
      dynamic_cast<LOCA::Bifurcation::HopfBord::ExtendedVector*>(results[i]);
    results_x[i] = &(hopfVecPtr->getXVec());
    results_y[i] = &(hopfVecPtr->getRealEigenVec());
    results_z[i] = &(hopfVecPtr->getImagEigenVec());
    results_w[i] = &(hopfVecPtr->getFrequency());
    results_params[i] = &(hopfVecPtr->getBifParam());

    tmp_real[i] = inputs_x[i]->clone(NOX::ShapeCopy);
    tmp_imag[i] = inputs_x[i]->clone(NOX::ShapeCopy);
  }

  // Set next to last components to deriv. w.r.t. parameter
  inputs_x[m] = derivResidualParamPtr;
  inputs_y[m] = derivRealEigenResidualParamPtr;
  inputs_z[m] = derivImagEigenResidualParamPtr;
  results_x[m] = inputs_x[m]->clone(NOX::ShapeCopy);
  results_y[m] = inputs_x[m]->clone(NOX::ShapeCopy);
  results_z[m] = inputs_x[m]->clone(NOX::ShapeCopy);
  tmp_real[m] = inputs_x[m]->clone(NOX::ShapeCopy);
  tmp_imag[m] = inputs_x[m]->clone(NOX::ShapeCopy);

  results_y[m+1] = inputs_y[m]->clone(NOX::ShapeCopy);
  results_z[m+1] = inputs_y[m]->clone(NOX::ShapeCopy);
  tmp_real[m+1] = inputs_x[m]->clone(NOX::ShapeCopy);
  tmp_imag[m+1] = inputs_x[m]->clone(NOX::ShapeCopy);

  // verify underlying Jacobian is valid
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Solve J*results_x = inputs_x
  finalStatus = grpPtr->applyJacobianInverseMulti(params, inputs_x, 
						  results_x, m+1);
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);
  

  // Compute tmp = (inputs_y + i*inputs_z) - (dJ(J+iwB)(y+iz)/dx)*results_x
  for (int i=0; i<m+1; i++) {
    status = grpPtr->computeDCeDxa(hopfXVec.getRealEigenVec(), 
				   hopfXVec.getImagEigenVec(),
				   hopfXVec.getFrequency(),
				   *results_x[i],
				   hopfFVec.getRealEigenVec(), 
				   hopfFVec.getImagEigenVec(),
				   *tmp_real[i],
				   *tmp_imag[i]);
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);

    tmp_real[i]->update(1.0, *inputs_y[i], -1.0);
    tmp_imag[i]->update(1.0, *inputs_z[i], -1.0);
  } 
  tmp_real[m+1] = minusMassTimesZPtr;
  tmp_imag[m+1] = massTimesYPtr;

  // verify underlying Jacobian is valid
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // verify underlying mass matrix is valid
  if (!grpPtr->isMassMatrix()) {
    status = grpPtr->computeMassMatrix();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Solve (J+iwB)*(results_y + i*results_z) = (tmp_real + i*tmp_imag)
  status = grpPtr->applyComplexInverseMulti(params, tmp_real, tmp_imag, 
					    hopfXVec.getFrequency(),
					    results_y, results_z, m+2);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  
  // Compute and set results
  double ltc = lTransNorm(*results_y[m+1]);
  double ltd = lTransNorm(*results_z[m+1]);
  double ltg = lTransNorm(*results_y[m]);
  double lth = lTransNorm(*results_z[m]);
  double denom_param = ltd*ltg - ltc*lth;
  double lte, ltf;
 
  for (int i=0; i<m; i++) {
    lte = lTransNorm(*results_y[i]);
    ltf = lTransNorm(*results_z[i]);
    *results_params[i] = 
      (ltc*(inputs_params[i] - ltf) + ltd*(lte - inputs_w[i])) / denom_param;
    *results_w[i] = (inputs_params[i] - ltf + lth*(*results_params[i]))/-ltd;

    results_x[i]->update(-*results_params[i], *results_x[m], 1.0);
    results_y[i]->update(-*results_params[i], *results_y[m], 
			 -*results_w[i], *results_y[m+1], 1.0);
    results_z[i]->update(-*results_params[i], *results_z[m], 
			 -*results_w[i], *results_z[m+1], 1.0);

    delete tmp_real[i];
    delete tmp_imag[i];
  }

  delete results_x[m];
  delete results_y[m];
  delete results_z[m];
  delete tmp_real[m];
  delete tmp_imag[m];

  delete results_y[m+1];
  delete results_z[m+1];

  delete [] inputs_x;
  delete [] inputs_y;
  delete [] inputs_z;
  delete [] inputs_w;
  delete [] inputs_params;
  delete [] results_x;
  delete [] results_y;
  delete [] results_z;
  delete [] results_w;
  delete [] results_params;
  delete [] tmp_real;
  delete [] tmp_imag;

  return finalStatus;
}

bool
LOCA::Bifurcation::HopfBord::ExtendedGroup::isF() const 
{
  return isValidF;
}

bool
LOCA::Bifurcation::HopfBord::ExtendedGroup::isJacobian() const 
{
  return isValidJacobian;
}

bool
LOCA::Bifurcation::HopfBord::ExtendedGroup::isGradient() const 
{
  return false;
}

bool
LOCA::Bifurcation::HopfBord::ExtendedGroup::isNewton() const 
{
  return isValidNewton;
}
  
const NOX::Abstract::Vector&
LOCA::Bifurcation::HopfBord::ExtendedGroup::getX() const 
{
  return hopfXVec;
}

const NOX::Abstract::Vector&
LOCA::Bifurcation::HopfBord::ExtendedGroup::getF() const 
{
  return hopfFVec;
}

double
LOCA::Bifurcation::HopfBord::ExtendedGroup::getNormF() const 
{
  return hopfFVec.norm();
}

const NOX::Abstract::Vector&
LOCA::Bifurcation::HopfBord::ExtendedGroup::getGradient() const 
{
  LOCA::ErrorCheck::throwError(
		 "LOCA::Bifurcation::HopfBord::ExtendedGroup::getGradient()",
		 " - not implemented");
  return getNewton();
}

const NOX::Abstract::Vector&
LOCA::Bifurcation::HopfBord::ExtendedGroup::getNewton() const 
{
  return hopfNewtonVec;
}

double
LOCA::Bifurcation::HopfBord::ExtendedGroup::getNormNewtonSolveResidual() const 
{
  string callingFunction = 
    "LOCA::Bifurcation::TPBord::ExtendedGroup::getNormNewtonSolveResidual()";
  NOX::Abstract::Group::ReturnType finalStatus;
  LOCA::Bifurcation::HopfBord::ExtendedVector residual = hopfFVec;
  
  finalStatus = applyJacobian(hopfNewtonVec,residual);
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  residual.update(1.0, hopfFVec, 1.0);
  return residual.norm();
}

void
LOCA::Bifurcation::HopfBord::ExtendedGroup::printSolution(const double conParam) const 
{
  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "LOCA::Bifurcation::HopfBord::ExtendedGroup::printSolution\n";

    cout << "\tPrinting Solution Vector for conParam = " 
	 << LOCA::Utils::sci(conParam) << endl;
  }
  grpPtr->printSolution(conParam);

  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "\tPrinting Real Component of Eigenvector for bif param = " 
	 << LOCA::Utils::sci(getBifParam()) << endl;
  }

  grpPtr->printSolution(hopfXVec.getRealEigenVec(), hopfXVec.getBifParam());
  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "\tPrinting Imaginary Component of Eigenvector for frequency = " 
	 << LOCA::Utils::sci(hopfXVec.getFrequency()) << endl;
  }
  grpPtr->printSolution(hopfXVec.getImagEigenVec(), 
			hopfXVec.getFrequency());
}

void
LOCA::Bifurcation::HopfBord::ExtendedGroup::printSolution(
					     const NOX::Abstract::Vector& x_,
					     const double conParam) const 
{
  const LOCA::Bifurcation::HopfBord::ExtendedVector& hopf_x = 
    dynamic_cast<const LOCA::Bifurcation::HopfBord::ExtendedVector&>(x_);
  
  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "LOCA::Bifurcation::HopfBord::ExtendedGroup::printSolution\n";

    cout << "\tPrinting Solution Vector for conParam = " 
	 << LOCA::Utils::sci(conParam) << endl;
  }
  grpPtr->printSolution(hopf_x.getXVec(), conParam);

  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "\tPrinting Real Component of Eigenvector for bif param = " 
	 << LOCA::Utils::sci(hopf_x.getBifParam()) << endl;
  }

  grpPtr->printSolution(hopf_x.getRealEigenVec(), hopf_x.getBifParam());
  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "\tPrinting Imaginary Component of Eigenvector for frequency = " 
	 << LOCA::Utils::sci(hopf_x.getFrequency()) << endl;
  }
  grpPtr->printSolution(hopf_x.getImagEigenVec(), 
			hopf_x.getFrequency());
}

const LOCA::Continuation::AbstractGroup&
LOCA::Bifurcation::HopfBord::ExtendedGroup::getUnderlyingGroup() const
{
  return *grpPtr;
}

LOCA::Continuation::AbstractGroup&
LOCA::Bifurcation::HopfBord::ExtendedGroup::getUnderlyingGroup()
{
  return *grpPtr;
}

void
LOCA::Bifurcation::HopfBord::ExtendedGroup::init(bool perturbSoln, 
						 double perturbSize)
{
  double ldy;
  double ldz;
  double denom;
  double a;
  double b;

  hopfXVec.getBifParam() = getBifParam();

  // Rescale and rotate complex eigenvector so normalization condition is met
  ldy = lTransNorm(hopfXVec.getRealEigenVec());
  ldz = lTransNorm(hopfXVec.getImagEigenVec());

  if (ldy == 0.0) {
     LOCA::ErrorCheck::throwError(
       "LOCA::Bifurcation::HopfBord::ExtendedGroup::init()",
       "Real component of eigenvector cannot be orthogonal to length vector!");
  }

  denom = ldy*ldy + ldz*ldz;
  a =  ldy/denom;
  b = -ldz/denom;
  NOX::Abstract::Vector* y_tmp = hopfXVec.getRealEigenVec().clone();

  // y <- a*y - b*z
  hopfXVec.getRealEigenVec().update(-b, hopfXVec.getImagEigenVec(), a);

  // z <- a*z + b*y
  hopfXVec.getImagEigenVec().update(b, *y_tmp, a);

  delete y_tmp;

  if (perturbSoln) {
    if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "\tIn LOCA::Bifurcation::HopfBord::ExtendedGroup::init(), applying random perturbation to initial solution of size:" 
	 << LOCA::Utils::sci(perturbSize) << endl;
    }
    NOX::Abstract::Vector *perturb = hopfXVec.getXVec().clone(NOX::ShapeCopy);
    perturb->random();
    perturb->scale(hopfXVec.getXVec());
    hopfXVec.getXVec().update(perturbSize, *perturb, 1.0);
    grpPtr->setX(hopfXVec.getXVec());
    delete perturb;
  }
}

double
LOCA::Bifurcation::HopfBord::ExtendedGroup::lTransNorm(
					const NOX::Abstract::Vector& n) const
{
  return lengthVecPtr->innerProduct(n) / lengthVecPtr->length();
}
