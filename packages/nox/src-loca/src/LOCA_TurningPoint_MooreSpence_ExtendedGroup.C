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

#include "LOCA_TurningPoint_MooreSpence_ExtendedGroup.H"
#include "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"
#include "LOCA_TurningPoint_MooreSpence_SolverStrategy.H"
#include "LOCA_Parameter_Vector.H"
#include "NOX_Parameter_List.H"
#include "LOCA_GlobalData.H"
#include "LOCA_Factory.H"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_Utils.H"
#include "LOCA_ErrorCheck.H"

LOCA::TurningPoint::MooreSpence::ExtendedGroup::ExtendedGroup(
	 const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
         const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<NOX::Parameter::List>& tpParams,
	 const Teuchos::RefCountPtr<LOCA::TurningPoint::MooreSpence::AbstractGroup>& g)
  : LOCA::Extended::MultiAbstractGroup(),
    LOCA::MultiContinuation::AbstractGroup(),
    globalData(global_data),
    parsedParams(topParams),
    turningPointParams(tpParams),
    grpPtr(g),
    xMultiVec(g->getX(), 2),
    fMultiVec(g->getX(), 2),
    newtonMultiVec(g->getX(), 2),
    lengthMultiVec(),
    xVec(),
    fVec(),
    ffMultiVec(),
    dfdpMultiVec(),
    newtonVec(),
    solverStrategy(),
    index_f(1),
    index_dfdp(1),
    bifParamID(1), 
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false)
{
  const char *func = "LOCA::TurningPoint::MooreSpence::ExtendedGroup()";

  if (!turningPointParams->isParameter("Bifurcation Parameter")) {
    globalData->locaErrorCheck->throwError(func,
				 "\"Bifurcation Parameter\" name is not set!");
  }
  string bifParamName = turningPointParams->getParameter(
						  "Bifurcation Parameter",
						  "None");
  const ParameterVector& p = grpPtr->getParams();
  bifParamID[0] = p.getIndex(bifParamName);

  if (!turningPointParams->isParameter("Length Normalization Vector")) {
    globalData->locaErrorCheck->throwError(func,
			   "\"Length Normalization Vector\" is not set!");
  }
  Teuchos::RefCountPtr<NOX::Abstract::Vector> lenVecPtr = 
    (*turningPointParams).INVALID_TEMPLATE_QUALIFIER 
      getRcpParameter<NOX::Abstract::Vector>("Length Normalization Vector");

  if (!turningPointParams->isParameter("Initial Null Vector")) {
    globalData->locaErrorCheck->throwError(func,
				 "\"Initial Null Vector\" is not set!");
  }
  Teuchos::RefCountPtr<NOX::Abstract::Vector> nullVecPtr = 
    (*turningPointParams).INVALID_TEMPLATE_QUALIFIER 
      getRcpParameter<NOX::Abstract::Vector>("Initial Null Vector");

  bool perturbSoln = turningPointParams->getParameter(
					       "Perturb Initial Solution", 
					       false);
  double perturbSize = turningPointParams->getParameter(
						 "Relative Perturbation Size", 
						 1.0e-3);

  lengthMultiVec = 
    Teuchos::rcp(lenVecPtr->createMultiVector(1, NOX::DeepCopy));
  xMultiVec.getColumn(0).getNullVec() = *nullVecPtr;

  // Instantiate solver strategy
  solverStrategy = 
    globalData->locaFactory->createMooreSpenceSolverStrategy(
							  parsedParams,
							  turningPointParams);

  // Set up multi-vector views
  setupViews(); 

  init(perturbSoln, perturbSize);
}

LOCA::TurningPoint::MooreSpence::ExtendedGroup::ExtendedGroup(
		const LOCA::TurningPoint::MooreSpence::ExtendedGroup& source, 
		NOX::CopyType type)
  : globalData(source.globalData),
    parsedParams(source.parsedParams),
    turningPointParams(source.turningPointParams),
    grpPtr(Teuchos::rcp(dynamic_cast<LOCA::TurningPoint::MooreSpence::AbstractGroup*>(source.grpPtr->clone(type)))),
    xMultiVec(source.xMultiVec, type),
    fMultiVec(source.fMultiVec, type),
    newtonMultiVec(source.newtonMultiVec, type),
    lengthMultiVec(Teuchos::rcp(source.lengthMultiVec->clone(type))),
    xVec(),
    fVec(),
    ffMultiVec(),
    dfdpMultiVec(),
    newtonVec(),
    lengthVec(),
    solverStrategy(source.solverStrategy),
    index_f(1),
    index_dfdp(1),
    bifParamID(source.bifParamID),
    isValidF(source.isValidF),
    isValidJacobian(source.isValidJacobian),
    isValidNewton(source.isValidNewton) 
{

  // Instantiate solver strategy
  solverStrategy = 
    globalData->locaFactory->createMooreSpenceSolverStrategy(
						      parsedParams,
						      turningPointParams);

  // Set up multi-vector views
  setupViews();

  if (type == NOX::ShapeCopy) {
    isValidF = false;
    isValidJacobian = false;
    isValidNewton = false;
  }
}

LOCA::TurningPoint::MooreSpence::ExtendedGroup::~ExtendedGroup() 
{
}

LOCA::TurningPoint::MooreSpence::ExtendedGroup&
LOCA::TurningPoint::MooreSpence::ExtendedGroup::operator=(
		 const LOCA::TurningPoint::MooreSpence::ExtendedGroup& source) 
{

  // Protect against A = A
  if (this != &source) {
    
    // Copy values
    globalData = source.globalData;
    parsedParams = source.parsedParams;
    turningPointParams = source.turningPointParams;
    *grpPtr = *source.grpPtr;
    xMultiVec = source.xMultiVec;
    fMultiVec = source.fMultiVec;
    newtonMultiVec = source.newtonMultiVec;
    *lengthMultiVec = *source.lengthMultiVec;
    index_f = source.index_f;
    index_dfdp = source.index_dfdp;
    bifParamID = source.bifParamID;
    isValidF = source.isValidF;
    isValidJacobian = source.isValidJacobian;
    isValidNewton = source.isValidNewton;

    // set up views again just to be safe
    setupViews();

    // Instantiate solver strategy
    solverStrategy = 
      globalData->locaFactory->createMooreSpenceSolverStrategy(
				   parsedParams,
				   turningPointParams);
  }

  return *this;
}

NOX::Abstract::Group&
LOCA::TurningPoint::MooreSpence::ExtendedGroup::operator=(
					   const NOX::Abstract::Group& source)
{
  return *this = 
   dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedGroup&>(source);
}

NOX::Abstract::Group*
LOCA::TurningPoint::MooreSpence::ExtendedGroup::clone(
						    NOX::CopyType type) const 
{
  return new LOCA::TurningPoint::MooreSpence::ExtendedGroup(*this, type);
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setX(
					      const NOX::Abstract::Vector& y) 
{
  const LOCA::TurningPoint::MooreSpence::ExtendedVector& yy = 
    dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedVector&>(y);
  grpPtr->setX( yy.getXVec() );
  *xVec = y;
  setBifParam(xVec->getBifParam());

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeX(
					      const NOX::Abstract::Group& g, 
					      const NOX::Abstract::Vector& d,
					      double step) 
{
  const LOCA::TurningPoint::MooreSpence::ExtendedGroup& gg = 
    dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedGroup&>(g);
  const LOCA::TurningPoint::MooreSpence::ExtendedVector& dd = 
    dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedVector&>(d);

  grpPtr->computeX(*(gg.grpPtr), dd.getXVec(), step);
  xVec->update(1.0, gg.getX(), step, dd, 0.0);
  setBifParam(xVec->getBifParam());

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeF() 
{
  if (isValidF)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeF()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Compute underlying F
  if (!grpPtr->isF()) {
    status = grpPtr->computeF();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }
  fVec->getXVec() = grpPtr->getF();
  
  // Compute underlying Jacobian
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // Compute J*n
  status = grpPtr->applyJacobian(xVec->getNullVec(), 
				 fVec->getNullVec());
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  
  // Compute phi^T*n
  fVec->getBifParam() = lTransNorm(xVec->getNullVec()) - 1.0;
  
  isValidF = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeJacobian() 
{
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeJacobian()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Compute underlying df/dp (may invalidate underlying data)
  status = grpPtr->computeDfDpMulti(bifParamID, 
				    fMultiVec.getXMultiVec(), 
				    isValidF);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute underlying dJn/dp (may invalidate underlying data)
  status = grpPtr->computeDJnDpMulti(bifParamID,
				     xVec->getNullVec(), 
				     fMultiVec.getNullMultiVec(), 
				     isValidF);

  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  solverStrategy->setBlocks(
		  grpPtr, 
		  Teuchos::rcp(this, false),
		  Teuchos::rcp(&(xVec->getNullVec()), false), 
		  Teuchos::rcp(&(fVec->getNullVec()), false),
		  Teuchos::rcp(&(fMultiVec.getColumn(1).getXVec()), false), 
		  Teuchos::rcp(&(fMultiVec.getColumn(1).getNullVec()), false));

  isValidJacobian = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeGradient() 
{
  return NOX::Abstract::Group::NotDefined;
}
   
NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeNewton(
						 NOX::Parameter::List& params) 
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeNewton()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Make sure F is valid
  if (!isF()) {
    status = computeF();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }
  
  // Make sure Jacobian is valid
  if (!isJacobian()) {
    status = computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // zero out newton vec -- used as initial guess for some linear solvers
  newtonMultiVec.init(0.0);

  // solve using contiguous
  status = solverStrategy->solve(params, fMultiVec, newtonMultiVec, true);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  newtonVec->scale(-1.0);

  isValidNewton = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobian(
					  const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  // Convert input, result to multivectors
  NOX::Abstract::MultiVector* mv_input = 
    input.createMultiVector(1, NOX::DeepCopy);
  NOX::Abstract::MultiVector* mv_result = 
    result.createMultiVector(1, NOX::DeepCopy);

  // Call multivector version of applyJacobian
  NOX::Abstract::Group::ReturnType status = 
    applyJacobianMultiVector(*mv_input, *mv_result);

  // Copy result
  result = (*mv_result)[0];

  // delete temporary multivectors
  delete mv_input;
  delete mv_result;

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianTranspose(
					  const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  // Convert input, result to multivectors
  NOX::Abstract::MultiVector* mv_input = 
    input.createMultiVector(1, NOX::DeepCopy);
  NOX::Abstract::MultiVector* mv_result = 
    result.createMultiVector(1, NOX::DeepCopy);

  // Call multivector version of applyJacobianTranspose
  NOX::Abstract::Group::ReturnType status = 
    applyJacobianTransposeMultiVector(*mv_input, *mv_result);

  // Copy result
  result = (*mv_result)[0];

  // delete temporary multivectors
  delete mv_input;
  delete mv_result;

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianInverse(
					  NOX::Parameter::List& params, 
					  const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  // Convert input, result to multivectors
  NOX::Abstract::MultiVector* mv_input = 
    input.createMultiVector(1, NOX::DeepCopy);
  NOX::Abstract::MultiVector* mv_result = 
    result.createMultiVector(1, NOX::DeepCopy);

  // Call multivector version of applyJacobianInverse
  NOX::Abstract::Group::ReturnType status = 
    applyJacobianInverseMultiVector(params, *mv_input, *mv_result);

  // Copy result
  result = (*mv_result)[0];

  // delete temporary multivectors
  delete mv_input;
  delete mv_result;

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianMultiVector(
				     const NOX::Abstract::MultiVector& input,
				     NOX::Abstract::MultiVector& result) const 
{
  string callingFunction = 
    "LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianMultiVector()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  if (!isJacobian()) {
    globalData->locaErrorCheck->throwError(callingFunction,
				 "Called with invalid Jacobian!");
  }

  // Cast vectors to TP vectors
  const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& tp_input = 
    dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector&>(input);
  LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& tp_result = 
    dynamic_cast<LOCA::TurningPoint::MooreSpence::ExtendedMultiVector&>(result);

  // Get constant references to input vector components
  const NOX::Abstract::MultiVector& input_x = tp_input.getXMultiVec();
  const NOX::Abstract::MultiVector& input_null = tp_input.getNullMultiVec();
  NOX::Abstract::MultiVector::DenseMatrix input_param = 
    tp_input.getScalars();

  // Get non-constant references to result vector components
  NOX::Abstract::MultiVector& result_x = tp_result.getXMultiVec();
  NOX::Abstract::MultiVector& result_null = tp_result.getNullMultiVec();
  NOX::Abstract::MultiVector::DenseMatrix result_param = 
    tp_result.getScalars();

  // Temporary vector
  NOX::Abstract::MultiVector *tmp = input_null.clone(NOX::ShapeCopy);

  // verify underlying Jacobian is valid
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  // compute J*x
  status = grpPtr->applyJacobianMultiVector(input_x, result_x);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // compute J*x + p*dR/dp
  result_x.update(Teuchos::NO_TRANS, 1.0, dfdpMultiVec->getXMultiVec(), 
		  input_param);

  // compute J*y
  status = grpPtr->applyJacobianMultiVector(input_null, result_null);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // compute J*y + p*dJy/dp
  result_null.update(Teuchos::NO_TRANS, 1.0, dfdpMultiVec->getNullMultiVec(), 
		     input_param);

  // compute (dJy/dx)*x
  status = grpPtr->computeDJnDxaMulti(xVec->getNullVec(), fVec->getNullVec(),
				      input_x, *tmp);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  
  // compute (dJy/dx)*x + J*y + p*dJy/dp
  result_null.update(1.0, *tmp, 1.0);

  // compute l^T*y
  lTransNorm(input_null, result_param);

  delete tmp;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianTransposeMultiVector(
				     const NOX::Abstract::MultiVector& input,
				     NOX::Abstract::MultiVector& result) const 
{
  globalData->locaErrorCheck->throwError(
		  "LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianTransposeMultiVector()",
		  "Method not implemented!");

  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianInverseMultiVector( 
				     NOX::Parameter::List& params, 	
				     const NOX::Abstract::MultiVector& input,
				     NOX::Abstract::MultiVector& result) const 
{
  // Cast vectors to TP vectors
  const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& tp_input = 
    dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector&>(input);
  LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& tp_result = 
    dynamic_cast<LOCA::TurningPoint::MooreSpence::ExtendedMultiVector&>(result);
  
  return solverStrategy->solve(params, tp_input, tp_result, false);
}

bool
LOCA::TurningPoint::MooreSpence::ExtendedGroup::isF() const 
{
  return isValidF;
}

bool
LOCA::TurningPoint::MooreSpence::ExtendedGroup::isJacobian() const 
{
  return isValidJacobian;
}

bool
LOCA::TurningPoint::MooreSpence::ExtendedGroup::isGradient() const 
{
  return false;
}

bool
LOCA::TurningPoint::MooreSpence::ExtendedGroup::isNewton() const 
{
  return isValidNewton;
}
  
const NOX::Abstract::Vector&
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getX() const 
{
  return *xVec;
}

const NOX::Abstract::Vector&
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getF() const 
{
  return *fVec;
}

double
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getNormF() const 
{
  return fVec->norm();
}

const NOX::Abstract::Vector&
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getGradient() const 
{
  globalData->locaErrorCheck->throwError(
	      "LOCA::TurningPoint::MooreSpence::ExtendedGroup::getGradient()",
	      " - not implemented");
  return getNewton();
}

const NOX::Abstract::Vector&
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getNewton() const 
{
  return *newtonVec;
}

double
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getNormNewtonSolveResidual() const 
{
  string callingFunction = 
    "LOCA::TurningPoint::MooreSpence::ExtendedGroup::getNormNewtonSolveResidual()";
  NOX::Abstract::Group::ReturnType finalStatus;
  LOCA::TurningPoint::MooreSpence::ExtendedVector residual = *fVec;
  
  finalStatus = applyJacobian(*newtonVec, residual);
  globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

  residual.update(1.0, *fVec, 1.0);
  return residual.norm();
}

LOCA::Extended::MultiAbstractGroup&
LOCA::TurningPoint::MooreSpence::ExtendedGroup::operator=(
			const LOCA::Extended::MultiAbstractGroup& source)
{
  return *this = 
   dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedGroup&>(source);
}

Teuchos::RefCountPtr<const LOCA::MultiContinuation::AbstractGroup>
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getUnderlyingGroup() const
{
  return grpPtr;
}

Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getUnderlyingGroup()
{
  return grpPtr;
}

LOCA::MultiContinuation::AbstractGroup&
LOCA::TurningPoint::MooreSpence::ExtendedGroup::operator=(
			const LOCA::MultiContinuation::AbstractGroup& source)
{
  return *this = 
   dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedGroup&>(source);
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setParamsMulti(
			  const vector<int>& paramIDs, 
			  const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  grpPtr->setParamsMulti(paramIDs, vals);
  for (unsigned int i=0; i<paramIDs.size(); i++) {
    if (paramIDs[i] == bifParamID[0])
      setBifParam(vals(i,0));
  }
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeDfDpMulti(
					    const vector<int>& paramIDs, 
					    NOX::Abstract::MultiVector& dfdp, 
					    bool isValid_F)
{
   string callingFunction = 
    "LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeDfDpMulti()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Cast dfdp to TP vector
  LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& tp_dfdp = 
    dynamic_cast<LOCA::TurningPoint::MooreSpence::ExtendedMultiVector&>(dfdp);

  // Compute df/dp
  status = grpPtr->computeDfDpMulti(paramIDs, tp_dfdp.getXMultiVec(),
				    isValid_F);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute d(Jn)/dp
  status = grpPtr->computeDJnDpMulti(paramIDs, xVec->getNullVec(),
				     tp_dfdp.getNullMultiVec(), isValid_F);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Set parameter components
  if (!isValid_F)
    tp_dfdp.getScalar(0,0) = lTransNorm(xVec->getNullVec());
  for (int i=0; i<dfdp.numVectors()-1; i++)
    tp_dfdp.getScalar(0,i+1) = 0.0;

  return finalStatus;
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::projectToDraw(
					       const NOX::Abstract::Vector& x,
					       double *px) const
{
  const LOCA::TurningPoint::MooreSpence::ExtendedVector& mx = 
    dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedVector&>(x);

  grpPtr->projectToDraw(mx.getXVec(), px);
  px[grpPtr->projectToDrawDimension()] = mx.getBifParam();
}

int
LOCA::TurningPoint::MooreSpence::ExtendedGroup::projectToDrawDimension() const
{
  return grpPtr->projectToDrawDimension() + 1;
}

LOCA::Continuation::AbstractGroup&
LOCA::TurningPoint::MooreSpence::ExtendedGroup::operator=(
			const LOCA::Continuation::AbstractGroup& source)
{
  return *this = 
   dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedGroup&>(source);
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setParams(
					      const LOCA::ParameterVector& p) 
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  grpPtr->setParams(p);
  setBifParam(p[bifParamID[0]]);
}

const LOCA::ParameterVector&
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getParams() const 
{
  return grpPtr->getParams();
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setParam(int paramID, 
							 double val)
{
  if (paramID == bifParamID[0])
    setBifParam(val);
  else
    grpPtr->setParam(paramID, val);
}

double
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getParam(int paramID) const
{
  return grpPtr->getParam(paramID);
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setParam(string paramID, 
							 double val)
{
  const LOCA::ParameterVector& pVec = grpPtr->getParams();
  if (pVec.getIndex(paramID) == bifParamID[0])
    setBifParam(val);
  else
    grpPtr->setParam(paramID, val);
}

double
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getParam(string paramID) const
{
  return grpPtr->getParam(paramID);
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeDfDp(int paramID, 
					      NOX::Abstract::Vector& result)
{
  vector<int> paramIDs(1);
  paramIDs[0] = paramID;
  NOX::Abstract::MultiVector *dfdp = result.createMultiVector(2);
  if (isValidF)
    (*dfdp)[0] = *fVec;
  NOX::Abstract::Group::ReturnType status = computeDfDpMulti(paramIDs, *dfdp,
							     isValidF);
  result = (*dfdp)[1];
  delete dfdp;

  return status;
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::printSolution(
						  const double conParam) const 
{
  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "LOCA::TurningPoint::MooreSpence::ExtendedGroup::printSolution\n";

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
  grpPtr->printSolution(xVec->getNullVec(), xVec->getBifParam());
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::printSolution(
					     const NOX::Abstract::Vector& x_,
					     const double conParam) const 
{
  const LOCA::TurningPoint::MooreSpence::ExtendedVector& tp_x = 
    dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedVector&>(x_);

  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "LOCA::TurningPoint::MooreSpence::ExtendedGroup::printSolution\n";

    cout << "Turning Point located at: " << LOCA::Utils::sci(conParam) 
	 << "   " 
	 << LOCA::Utils::sci(tp_x.getBifParam())
	 << endl;

    cout << "\tPrinting Solution Vector for conParam = " 
	 << LOCA::Utils::sci(conParam) << endl;
  }
  grpPtr->printSolution(tp_x.getXVec(), conParam);
  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "\tPrinting Null Vector for bif param = " 
	 << LOCA::Utils::sci(tp_x.getBifParam()) << endl;
  }
  grpPtr->printSolution(tp_x.getNullVec(), tp_x.getBifParam());
}

double
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getBifParam() const 
{
  return grpPtr->getParam(bifParamID[0]);
}

double
LOCA::TurningPoint::MooreSpence::ExtendedGroup::lTransNorm(
					const NOX::Abstract::Vector& n) const
{
  return lengthVec->dot(n) / lengthVec->length();
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::lTransNorm(
			const NOX::Abstract::MultiVector& n,
			NOX::Abstract::MultiVector::DenseMatrix& result) const
{
  n.multiply(1.0 / lengthVec->length(), *lengthMultiVec, result);
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setBifParam(double param) 
{
  grpPtr->setParam(bifParamID[0], param);
  xVec->getBifParam() = param;

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setupViews()
{
  index_f[0] = 0;
  index_dfdp[0] = 1;
  
  xVec = Teuchos::rcp(dynamic_cast<LOCA::TurningPoint::MooreSpence::ExtendedVector*>(&xMultiVec[0]), false);
  fVec = Teuchos::rcp(dynamic_cast<LOCA::TurningPoint::MooreSpence::ExtendedVector*>(&fMultiVec[0]), false);
  newtonVec = Teuchos::rcp(dynamic_cast<LOCA::TurningPoint::MooreSpence::ExtendedVector*>(&newtonMultiVec[0]), false);
  lengthVec = Teuchos::rcp(&(*lengthMultiVec)[0], false);

  ffMultiVec = Teuchos::rcp(dynamic_cast<LOCA::TurningPoint::MooreSpence::ExtendedMultiVector*>(fMultiVec.subView(index_f)));

  dfdpMultiVec = Teuchos::rcp(dynamic_cast<LOCA::TurningPoint::MooreSpence::ExtendedMultiVector*>(fMultiVec.subView(index_dfdp)));

}

void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::init(bool perturbSoln, 
						     double perturbSize)
{
  xVec->getBifParam() = getBifParam();

  // Rescale length vector so that the normalization condition is met
  double lVecDotNullVec = lTransNorm(xVec->getNullVec());

  if (lVecDotNullVec == 0.0) {
    globalData->locaErrorCheck->throwError(
		   "LOCA::TurningPoint::MooreSpence::ExtendedGroup::init()",
		   "null vector can be orthogonal to length-scaling vector");
  }
  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "\tIn LOCA::TurningPoint::MooreSpence::ExtendedGroup::init(), scaling null vector by:" 
	 << LOCA::Utils::sci(1.0 / lVecDotNullVec) << endl;
  }
  xVec->getNullVec().scale(1.0/lVecDotNullVec);

  if (perturbSoln) {
    if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "\tIn LOCA::TurningPoint::MooreSpence::ExtendedGroup::init(), applying random perturbation to initial solution of size:" 
	 << LOCA::Utils::sci(perturbSize) << endl;
    }
    NOX::Abstract::Vector *perturb = xVec->getXVec().clone(NOX::ShapeCopy);
    perturb->random();
    perturb->scale(xVec->getXVec());
    xVec->getXVec().update(perturbSize, *perturb, 1.0);
    grpPtr->setX(xVec->getXVec());
    delete perturb;
  }
}

