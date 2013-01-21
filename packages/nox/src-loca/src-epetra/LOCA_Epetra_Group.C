// $Id: LOCA_Epetra_Group.C,v 1.53 2009/03/25 18:50:28 etphipp Exp $ 
// $Source: /space/CVS/Trilinos/packages/nox/src-loca/src-epetra/LOCA_Epetra_Group.C,v $ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "LOCA_Epetra_Group.H"	          // class definition

#include "LOCA_Epetra_Interface_Required.H"        // class data members
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "LOCA_Epetra_ShiftInvertOperator.H"
#include "AztecOO.h"
#include "NOX_Epetra_LinearSystem_AztecOO.H"
#include "LOCA_Epetra_TransposeLinearSystem_AbstractStrategy.H"
#include "LOCA_Epetra_TransposeLinearSystem_Factory.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "NOX_Epetra_MultiVector.H"
#include "LOCA_Parameter_SublistParser.H"

#ifdef HAVE_NOX_EPETRAEXT
#include "EpetraExt_BlockCrsMatrix.h"
#include "EpetraExt_BlockVector.h"
#include "EpetraExt_BlockMultiVector.h"
#endif

LOCA::Epetra::Group::Group(
	    const Teuchos::RCP<LOCA::GlobalData>& global_data,
	    Teuchos::ParameterList& printingParams, 
	    const Teuchos::RCP<LOCA::Epetra::Interface::Required>& i, 
	    NOX::Epetra::Vector& initialGuess,
	    const LOCA::ParameterVector& p) :
  NOX::Epetra::Group(printingParams, i, initialGuess),
  LOCA::Abstract::Group(global_data),
  globalData(global_data),
  printParams(printingParams),
  params(p),
  userInterface(i),
  userInterfaceTime(),
  userInterfaceTimeMF(),
  userInterfaceFreeEnergy(),
  shiftedSharedLinearSystem(),
  isValidShiftedPrec(false),
  alpha_(1.0),
  beta_(0.0),
  tmpVectorPtr2(),
  scaleVecPtr(),
  tls_strategy(),
  complexSharedLinearSystem(),
  complexMatrix(),
  complexVec(),
  isValidComplex(false),
  isValidComplexPrec(false),
  separateMatrixMemoryDeclared(false)
{
}

LOCA::Epetra::Group::Group(
	    const Teuchos::RCP<LOCA::GlobalData>& global_data,
	    Teuchos::ParameterList& printingParams, 
	    const Teuchos::RCP<LOCA::Epetra::Interface::Required>& i, 
	    NOX::Epetra::Vector& initialGuess, 
	    const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys,
	    const LOCA::ParameterVector& p) :
  NOX::Epetra::Group(printingParams, i, initialGuess, linSys),
  LOCA::Abstract::Group(global_data),
  globalData(global_data),
  printParams(printingParams),
  params(p),
  userInterface(i),
  userInterfaceTime(),
  userInterfaceTimeMF(),
  userInterfaceFreeEnergy(),
  shiftedSharedLinearSystem(),
  isValidShiftedPrec(false),
  alpha_(1.0),
  beta_(0.0),
  tmpVectorPtr2(),
  scaleVecPtr(),
  tls_strategy(),
  complexSharedLinearSystem(),
  complexMatrix(),
  complexVec(),
  isValidComplex(false),
  isValidComplexPrec(false),
  separateMatrixMemoryDeclared(false)
{
}

LOCA::Epetra::Group::Group(
	const Teuchos::RCP<LOCA::GlobalData>& global_data,
	Teuchos::ParameterList& printingParams, 
	const Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent>& i, 
	NOX::Epetra::Vector& initialGuess, 
	const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys,
	const Teuchos::RCP<NOX::Epetra::LinearSystem>& shiftedLinSys,
	const LOCA::ParameterVector& p) :
  NOX::Epetra::Group(printingParams, i, initialGuess, linSys),
  LOCA::Abstract::Group(global_data),
  globalData(global_data),
  printParams(printingParams),
  params(p),
  userInterface(i), 
  userInterfaceTime(i),
  userInterfaceTimeMF(),
  userInterfaceFreeEnergy(),
  shiftedSharedLinearSystem(),
  isValidShiftedPrec(false),
  alpha_(1.0),
  beta_(0.0),
  tmpVectorPtr2(),
  scaleVecPtr(),
  tls_strategy(),
  complexSharedLinearSystem(),
  complexMatrix(),
  complexVec(),
  isValidComplex(false),
  isValidComplexPrec(false),
  separateMatrixMemoryDeclared(false)
{
  shiftedSharedLinearSystem = 
    Teuchos::rcp(new NOX::SharedObject<NOX::Epetra::LinearSystem, 
		                       LOCA::Epetra::Group>(shiftedLinSys));
}

LOCA::Epetra::Group::Group(
	    const Teuchos::RCP<LOCA::GlobalData>& global_data,
	    Teuchos::ParameterList& printingParams, 
	    const Teuchos::RCP<LOCA::Epetra::Interface::TimeDependentMatrixFree>& i, 
	    NOX::Epetra::Vector& initialGuess, 
	    const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys,
	    const Teuchos::RCP<NOX::Epetra::LinearSystem>& shiftedLinSys,
	    const LOCA::ParameterVector& p) :
  NOX::Epetra::Group(printingParams, i, initialGuess, linSys),
  LOCA::Abstract::Group(global_data),
  globalData(global_data),
  printParams(printingParams),
  params(p),
  userInterface(i),
  userInterfaceTime(),
  userInterfaceTimeMF(i),
  userInterfaceFreeEnergy(),
  shiftedSharedLinearSystem(),
  isValidShiftedPrec(false),
  alpha_(1.0),
  beta_(0.0),
  tmpVectorPtr2(),
  scaleVecPtr(),
  tls_strategy(),
  complexSharedLinearSystem(),
  complexMatrix(),
  complexVec(),
  isValidComplex(false),
  isValidComplexPrec(false),
  separateMatrixMemoryDeclared(false)
{
 shiftedSharedLinearSystem = 
    Teuchos::rcp(new NOX::SharedObject<NOX::Epetra::LinearSystem, 
		                       LOCA::Epetra::Group>(shiftedLinSys)); 
}

LOCA::Epetra::Group::Group(const LOCA::Epetra::Group& source, 
			   NOX::CopyType type) :
  NOX::Epetra::Group(source, type),
  LOCA::Abstract::Group(source, type),
  globalData(source.globalData),
  printParams(source.printParams),
  params(source.params),
  userInterface(source.userInterface),
  userInterfaceTime(source.userInterfaceTime),
  userInterfaceTimeMF(source.userInterfaceTimeMF),
  userInterfaceFreeEnergy(source.userInterfaceFreeEnergy),
  shiftedSharedLinearSystem(source.shiftedSharedLinearSystem),
  isValidShiftedPrec(source.isValidShiftedPrec),
  alpha_(source.alpha_),
  beta_(source.beta_),
  tmpVectorPtr2(),
  scaleVecPtr(),
  tls_strategy(),
  complexSharedLinearSystem(source.complexSharedLinearSystem),
  complexMatrix(source.complexMatrix),
  complexVec(source.complexVec),
  isValidComplex(source.isValidComplex),
  isValidComplexPrec(source.isValidComplexPrec),
  separateMatrixMemoryDeclared(source.separateMatrixMemoryDeclared)
{
  if (source.scaleVecPtr != Teuchos::null)
    scaleVecPtr = source.scaleVecPtr->clone(NOX::DeepCopy);

  // Take ownership of shared shifted system
  if (shiftedSharedLinearSystem != Teuchos::null && type == NOX::DeepCopy)
    shiftedSharedLinearSystem->getObject(this);

  // Take ownership of complex system
  if (complexSharedLinearSystem != Teuchos::null && type == NOX::DeepCopy)
    complexSharedLinearSystem->getObject(this);
}

LOCA::Epetra::Group::~Group() 
{
}

LOCA::Epetra::Group& 
LOCA::Epetra::Group::operator=(const LOCA::Epetra::Group& source)
{
  if (this != &source) {
    params = source.params;
    NOX::Epetra::Group::operator=(source);
    LOCA::Abstract::Group::copy(source);
    params = source.params;
    userInterface = source.userInterface;
    userInterfaceTime = source.userInterfaceTime;
    userInterfaceFreeEnergy = source.userInterfaceFreeEnergy;
    if (source.scaleVecPtr != Teuchos::null)
      scaleVecPtr = source.scaleVecPtr->clone(NOX::DeepCopy);
    
    // Take ownership of shared shifted system
    if (shiftedSharedLinearSystem != Teuchos::null)
      shiftedSharedLinearSystem->getObject(this);
    isValidShiftedPrec = source.isValidShiftedPrec;
    alpha_ = source.alpha_;
    beta_ = source.beta_;

    // Take ownership of complex system
    if (complexSharedLinearSystem != Teuchos::null)
      complexSharedLinearSystem->getObject(this);
    isValidComplex = source.isValidComplex;
    isValidComplexPrec = source.isValidComplexPrec;
    separateMatrixMemoryDeclared = source.separateMatrixMemoryDeclared;
  }
  return *this;
}

NOX::Abstract::Group& 
LOCA::Epetra::Group::operator=(const NOX::Abstract::Group& source)
{
  operator=(dynamic_cast<const Group&> (source));
  return *this;
}

NOX::Abstract::Group& 
LOCA::Epetra::Group::operator=(const NOX::Epetra::Group& source)
{
  operator=(dynamic_cast<const Group&> (source));
  return *this;
}

Teuchos::RCP<NOX::Abstract::Group>
LOCA::Epetra::Group::clone(NOX::CopyType type) const 
{
  return Teuchos::rcp(new LOCA::Epetra::Group(*this, type));
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::computeF() 
{

  if (isF())
    return Abstract::Group::Ok;
  
  // Set the parameters prior to computing F
  userInterface->setParameters(params);
  
  return NOX::Epetra::Group::computeF();
}

NOX::Abstract::Group::ReturnType 
LOCA::Epetra::Group::computeJacobian() 
{

  if (isJacobian())
    return Abstract::Group::Ok;
  
  // Set the parameters prior to computing F
  userInterface->setParameters(params);

  return NOX::Epetra::Group::computeJacobian();
}

NOX::Abstract::Group::ReturnType 
LOCA::Epetra::Group::applyJacobianTransposeInverse(
				    Teuchos::ParameterList& params, 
				    const NOX::Abstract::Vector& input, 
				    NOX::Abstract::Vector& result) const
{
  std::string callingFunction = 
    "LOCA::Epetra::Group::applyJacobianTransposeInverse()";
  NOX::Abstract::Group::ReturnType status;
  
  // Get non-const linsys
  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys = 
    sharedLinearSystem.getObject(this);

  // Get Jacobian operator
  Teuchos::RCP<Epetra_Operator> jac =
    linSys->getJacobianOperator();

  // Instantiate transpose solver
  LOCA::Epetra::TransposeLinearSystem::Factory tls_factory(globalData);
  if (tls_strategy == Teuchos::null)
    tls_strategy = tls_factory.create(Teuchos::rcp(&params, false), linSys);
     
  // Compute Jacobian transpose J^T
  tls_strategy->createJacobianTranspose();

  // Now compute preconditioner for J^T
  tls_strategy->createTransposePreconditioner(xVector, params);

  // Solve
  bool stat = 
    tls_strategy->applyJacobianTransposeInverse(
			  params, 
			  dynamic_cast<const NOX::Epetra::Vector&>(input),
			  dynamic_cast<NOX::Epetra::Vector&>(result));
  if (stat == true)
    status = NOX::Abstract::Group::Ok;
  else
    status = NOX::Abstract::Group::NotConverged;

  // Set original operators in linear system
  jac->SetUseTranspose(false);
  linSys->setJacobianOperatorForSolve(jac);
  linSys->destroyPreconditioner();

  return status;
}

NOX::Abstract::Group::ReturnType 
LOCA::Epetra::Group::applyJacobianTransposeInverseMultiVector(
				    Teuchos::ParameterList& params, 
				    const NOX::Abstract::MultiVector& input, 
				    NOX::Abstract::MultiVector& result) const
{
  std::string callingFunction = 
    "LOCA::Epetra::Group::applyJacobianTransposeInverseMultiVector()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // Get non-const linsys
  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys = 
    sharedLinearSystem.getObject(this);

  // Get Jacobian operator
  Teuchos::RCP<Epetra_Operator> jac =
    linSys->getJacobianOperator();

  // Instantiate transpose solver
  LOCA::Epetra::TransposeLinearSystem::Factory tls_factory(globalData);
  if (tls_strategy == Teuchos::null)
    tls_strategy = tls_factory.create(Teuchos::rcp(&params, false), linSys);
     
  // Compute Jacobian transpose J^T
  tls_strategy->createJacobianTranspose();

  // Now compute preconditioner for J^T
  tls_strategy->createTransposePreconditioner(xVector, params);

  // Solve for each RHS
  int m = input.numVectors();
  for (int i=0; i<m; i++) {
    bool stat = 
      tls_strategy->applyJacobianTransposeInverse(
			  params, 
			  dynamic_cast<const NOX::Epetra::Vector&>(input[i]),
			  dynamic_cast<NOX::Epetra::Vector&>(result[i]));
    if (stat == true)
      status = NOX::Abstract::Group::Ok;
    else
      status = NOX::Abstract::Group::NotConverged;
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  // Set original operators in linear system
  jac->SetUseTranspose(false);
  linSys->setJacobianOperatorForSolve(jac);
  linSys->destroyPreconditioner();

  return finalStatus;
}

void
LOCA::Epetra::Group::copy(const NOX::Abstract::Group& source)
{
  *this = source;
}

void 
LOCA::Epetra::Group::setParams(const LOCA::ParameterVector& p)
{
  resetIsValid();
  params = p;
}

void
LOCA::Epetra::Group::setParam(int paramID, double val)
{
  resetIsValid();
  params.setValue(paramID, val);
}

void
LOCA::Epetra::Group::setParam(std::string paramID, double val)
{
  resetIsValid();
  params.setValue(paramID, val);
}

const LOCA::ParameterVector& 
LOCA::Epetra::Group::getParams() const 
{
  return params;
}

double
LOCA::Epetra::Group::getParam(int paramID) const
{
  return params.getValue(paramID);
}

double
LOCA::Epetra::Group::getParam(std::string paramID) const
{
  return params.getValue(paramID);
}

void
LOCA::Epetra::Group:: preProcessContinuationStep(
			     LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  userInterface->preProcessContinuationStep(stepStatus, *this);
}

void
LOCA::Epetra::Group:: postProcessContinuationStep(
			     LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  userInterface->postProcessContinuationStep(stepStatus, *this);
}

void
LOCA::Epetra::Group::projectToDraw(const NOX::Abstract::Vector& x,
				   double *px) const
{
  const NOX::Epetra::Vector& ex = 
    dynamic_cast<const NOX::Epetra::Vector&>(x);
  userInterface->projectToDraw(ex, px);
}

int
LOCA::Epetra::Group::projectToDrawDimension() const
{
  return userInterface->projectToDrawDimension();
}

double
LOCA::Epetra::Group::computeScaledDotProduct(
				       const NOX::Abstract::Vector& a,
				       const NOX::Abstract::Vector& b) const
{
  if (scaleVecPtr == Teuchos::null)
    return a.innerProduct(b) / a.length();
  else {
    Teuchos::RCP<NOX::Abstract::Vector> as = a.clone(NOX::DeepCopy);
    Teuchos::RCP<NOX::Abstract::Vector> bs = b.clone(NOX::DeepCopy);
    double d;

    as->scale(*scaleVecPtr);
    bs->scale(*scaleVecPtr);
    d = as->innerProduct(*bs);

    return d;
  }
}

void
LOCA::Epetra::Group::printSolution(const double conParam) const
{
  printSolution(xVector, conParam);
}

void
LOCA::Epetra::Group::printSolution(const NOX::Abstract::Vector& x_,
				      const double conParam) const
{
  printSolution(dynamic_cast<const NOX::Epetra::Vector&>(x_), conParam);
}

void
LOCA::Epetra::Group::scaleVector(NOX::Abstract::Vector& x) const
{
  if (scaleVecPtr == Teuchos::null)
    x.scale(1.0 / sqrt(static_cast<double>(x.length())));
  else 
    x.scale(*scaleVecPtr);
}

NOX::Abstract::Group::ReturnType 
LOCA::Epetra::Group::augmentJacobianForHomotopy(double a, double b)
{

  //Allocate temporary vectors if not aready done
  if (tmpVectorPtr == Teuchos::null)
    tmpVectorPtr = Teuchos::rcp(new Epetra_Vector(xVector.getEpetraVector()));
  if (tmpVectorPtr2 == Teuchos::null)
    tmpVectorPtr2 = Teuchos::rcp(new Epetra_Vector(xVector.getEpetraVector()));

  tmpVectorPtr2->PutScalar(b);

  // See if it is an Epetra_CrsMatrix
  Teuchos::RCP<const Epetra_CrsMatrix> constTestCrs;
  Teuchos::RCP<Epetra_CrsMatrix> testCrs;
  constTestCrs = Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>
    (sharedLinearSystem.getObject(this)->getJacobianOperator());
  if (constTestCrs != Teuchos::null) {
    testCrs = Teuchos::rcp_const_cast<Epetra_CrsMatrix>(constTestCrs);
    testCrs->Scale(a);
    testCrs->ExtractDiagonalCopy(*tmpVectorPtr);
    tmpVectorPtr->Update(1.0, *tmpVectorPtr2, 1.0);
    testCrs->ReplaceDiagonalValues(*tmpVectorPtr);
    return LOCA::Abstract::Group::Ok;

  }

  // See if it is an Epetra_VbrMatrix
  Teuchos::RCP<const Epetra_VbrMatrix> constTestVbr;
  Teuchos::RCP<Epetra_VbrMatrix> testVbr;
  constTestVbr = Teuchos::rcp_dynamic_cast<const Epetra_VbrMatrix>
    (sharedLinearSystem.getObject(this)->getJacobianOperator());
  if (constTestVbr != Teuchos::null) {
    testVbr = Teuchos::rcp_const_cast<Epetra_VbrMatrix>(constTestVbr);
    testVbr->Scale(a);
    testVbr->ExtractDiagonalCopy(*tmpVectorPtr);
    tmpVectorPtr->Update(1.0, *tmpVectorPtr2, 1.0);
    testVbr->ReplaceDiagonalValues(*tmpVectorPtr);
    return LOCA::Abstract::Group::Ok;
  }

  // Otherwise this alg won't work -- return NotDefined

  return LOCA::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::computeShiftedMatrix(double alpha, double beta)
{
  // We store a real shifted matrix
  if (userInterfaceTime != Teuchos::null) {
    Teuchos::RCP<Epetra_Operator> mass = 
      shiftedSharedLinearSystem->getObject(this)->getJacobianOperator();

    bool res = 
      userInterfaceTime->computeShiftedMatrix(alpha, beta, 
					      xVector.getEpetraVector(), 
					      *mass);
    
    // Check if Jacobian and mass matrices are the same, in which case
    // the Jacobian is no longer valid
    Teuchos::RCP<Epetra_Operator> jac = 
      sharedLinearSystem.getObject(this)->getJacobianOperator();
    if (mass.get() == jac.get())
      isValidJacobian = false;

    if (res)
      return NOX::Abstract::Group::Ok;
    else
      return NOX::Abstract::Group::Failed;
  }

  // Matrix free shifted matrix
  else if (userInterfaceTimeMF != Teuchos::null) {
    alpha_ = alpha;
    beta_ = beta;
    return NOX::Abstract::Group::Ok;
  }

  // Failure
  else
    return NOX::Abstract::Group::BadDependency;
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyShiftedMatrix(const NOX::Abstract::Vector& input,
                                        NOX::Abstract::Vector& result) const
{
  const NOX::Epetra::Vector& epetra_input = 
    dynamic_cast<const NOX::Epetra::Vector&>(input);
  NOX::Epetra::Vector& epetra_result = 
    dynamic_cast<NOX::Epetra::Vector&>(result);

  // We store a real shifted matrix
  if (userInterfaceTime != Teuchos::null) {
    bool res =
      shiftedSharedLinearSystem->getObject()->applyJacobian(epetra_input, 
							    epetra_result);
    if (res)
      return NOX::Abstract::Group::Ok;
    else
      return NOX::Abstract::Group::Failed;
  }

  // For matrix-free mode, call interface for applying mass matrix
  else if (userInterfaceTimeMF != Teuchos::null) {
     bool res = userInterfaceTimeMF->applyShiftedMatrix(alpha_, beta_, 
							epetra_input,
							epetra_result);
     if (res)
      return NOX::Abstract::Group::Ok;
    else
      return NOX::Abstract::Group::Failed;
  }

  // Failure
  else
    return NOX::Abstract::Group::BadDependency;
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyShiftedMatrixMultiVector(
				     const NOX::Abstract::MultiVector& input,
				     NOX::Abstract::MultiVector& result) const
{
  std::string callingFunction = 
    "LOCA::Epetra::Group::applyShiftedMatrixMultiVector()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  for (int i=0; i<input.numVectors(); i++) {
    status = applyShiftedMatrix(input[i], result[i]);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyShiftedMatrixInverseMultiVector(
			        Teuchos::ParameterList& lsParams, 
				const NOX::Abstract::MultiVector& input,
				NOX::Abstract::MultiVector& result) const
{

  if (shiftedSharedLinearSystem != Teuchos::null) {

    NOX::Epetra::LinearSystem::PreconditionerReusePolicyType precPolicy = 
      sharedLinearSystem.getObject(this)->getPreconditionerPolicy();

    if (!isValidShiftedPrec) {

      if (precPolicy == NOX::Epetra::LinearSystem::PRPT_REBUILD) {
	shiftedSharedLinearSystem->getObject(this)->destroyPreconditioner();
	shiftedSharedLinearSystem->getObject(this)->
	  createPreconditioner(xVector, lsParams, false);
	isValidShiftedPrec = true;
      }
      else if (precPolicy == NOX::Epetra::LinearSystem::PRPT_RECOMPUTE) {
	sharedLinearSystem.getObject(this)->recomputePreconditioner(xVector, 
								    lsParams);
      }
      else if (precPolicy == NOX::Epetra::LinearSystem::PRPT_REUSE) {
	// Do Nothing!!!
      }

    }


    const NOX::Epetra::Vector* epetra_input;
    NOX::Epetra::Vector* epetra_result; 
    bool status;
    bool finalStatus = true;
    for (int i=0; i<input.numVectors(); i++) {
      epetra_input = dynamic_cast<const NOX::Epetra::Vector*>(&input[i]);
      epetra_result = dynamic_cast<NOX::Epetra::Vector*>(&result[i]);

      status = 
	shiftedSharedLinearSystem->getObject(this)->applyJacobianInverse(
							      lsParams, 
							      *epetra_input, 
							      *epetra_result);
      finalStatus = finalStatus && status;
    }
    
    if (finalStatus)
      return NOX::Abstract::Group::Ok;
    else
      return NOX::Abstract::Group::NotConverged;
  }

  else 
    return NOX::Abstract::Group::BadDependency;
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::computeSecondShiftedMatrix(double alpha, double beta)
{
  // Use memory from regular "jacobian" linear system for secondShiftedMatrix,
  // which is only valid if user declareSeparateMatrixMemory
  if (!separateMatrixMemoryDeclared) {
    globalData->locaErrorCheck->throwError(
       "LOCA::Epetra::Group::computeSecondShiftedMatrix()",
       "Group must be called with declareSeparateMatrixMemory() for use of second shifted matrix.");
    return NOX::Abstract::Group::BadDependency;
  }

  // We store a real shifted matrix
  if (userInterfaceTime != Teuchos::null) {
    Teuchos::RCP<Epetra_Operator> shift2 = 
      sharedLinearSystem.getObject(this)->getJacobianOperator();

    bool res = 
      userInterfaceTime->computeShiftedMatrix(alpha, beta, 
					      xVector.getEpetraVector(), 
					      *shift2);
    
    // Check if Jacobian and mass matrices are the same, in which case
    // the Jacobian is no longer valid
    isValidJacobian = false;

    if (res)
      return NOX::Abstract::Group::Ok;
    else
      return NOX::Abstract::Group::Failed;
  }

  // Failure
  else
    return NOX::Abstract::Group::BadDependency;
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applySecondShiftedMatrix(const NOX::Abstract::Vector& input,
                                        NOX::Abstract::Vector& result) const
{
  const NOX::Epetra::Vector& epetra_input = 
    dynamic_cast<const NOX::Epetra::Vector&>(input);
  NOX::Epetra::Vector& epetra_result = 
    dynamic_cast<NOX::Epetra::Vector&>(result);

  // We store a real shifted matrix
  if (userInterfaceTime != Teuchos::null) {
    bool res =
      sharedLinearSystem.getObject()->applyJacobian(epetra_input, 
						     epetra_result);
    if (res)
      return NOX::Abstract::Group::Ok;
    else
      return NOX::Abstract::Group::Failed;
  }

  // Failure
  else
    return NOX::Abstract::Group::BadDependency;
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applySecondShiftedMatrixMultiVector(
				     const NOX::Abstract::MultiVector& input,
				     NOX::Abstract::MultiVector& result) const
{
  std::string callingFunction = 
    "LOCA::Epetra::Group::applySecondShiftedMatrixMultiVector()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  for (int i=0; i<input.numVectors(); i++) {
    status = applySecondShiftedMatrix(input[i], result[i]);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  return finalStatus;
}

bool
LOCA::Epetra::Group::isComplex() const
{
  return isValidComplex;
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::computeComplex(double frequency)
{
  std::string callingFunction = 
    "LOCA::Epetra::Group::computeComplex()";

  // We must have the time-dependent interface
  if (userInterfaceTime == Teuchos::null)
    return NOX::Abstract::Group::BadDependency;

#ifdef HAVE_NOX_EPETRAEXT
  if (isValidComplex)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Get Jacobian matrix
  Teuchos::RCP<Epetra_RowMatrix> jac = Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(sharedLinearSystem.getObject(this)->getJacobianOperator());
   
  // Create complex matrix
  if (complexMatrix == Teuchos::null) {
    std::vector< std::vector<int> >rowStencil(2);
    std::vector<int> rowIndex;

    rowStencil[0].push_back(0); rowStencil[0].push_back(1);
    rowStencil[1].push_back(-1);  rowStencil[1].push_back(0);
    rowIndex.push_back(0); rowIndex.push_back(1);

    complexMatrix = Teuchos::rcp(new EpetraExt::BlockCrsMatrix(*jac,
							       rowStencil, 
							       rowIndex, 
							       jac->Comm()));

    // Construct global solution vector, the overlap vector, and importer between them
   complexVec = 
     Teuchos::rcp(new EpetraExt::BlockVector(jac->RowMatrixRowMap(), 
					     complexMatrix->RowMap()));  
  }

  // Get mass matrix M
  Teuchos::RCP<Epetra_RowMatrix> mass = Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(shiftedSharedLinearSystem->getObject(this)->getJacobianOperator());

  // Compute w*M
  status = computeShiftedMatrix(0.0, frequency);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							   finalStatus,
							   callingFunction);

  // Load w*M in complex matrix
  complexMatrix->LoadBlock(*mass, 1, 0);

  // Compute -w*M
  status = computeShiftedMatrix(0.0, -frequency);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							   finalStatus,
							   callingFunction);

  // Load -w*M in complex matrix
  complexMatrix->LoadBlock(*mass, 0, 1);

  // Compute J
  status = computeJacobian();
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							   finalStatus,
							   callingFunction);

  // Load J in complex matrix
  complexMatrix->LoadBlock(*jac, 0, 0);
  complexMatrix->LoadBlock(*jac, 1, 1);

  // Create linear system
  if (complexSharedLinearSystem == Teuchos::null) {

    NOX::Epetra::Vector nev(complexVec, NOX::Epetra::Vector::CreateView);

    Teuchos::RCP<Teuchos::ParameterList> lsParams = 
      globalData->parsedParams->getSublist("Linear Solver");

    // Create the Linear System
    Teuchos::RCP<NOX::Epetra::Interface::Required> iReq;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac;
    Teuchos::RCP<NOX::Epetra::LinearSystem> complexLinSys = 
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, 
							*lsParams,
							iReq, 
							iJac, 
							complexMatrix, 
							nev));
    complexLinSys->setJacobianOperatorForSolve(complexMatrix);
    complexSharedLinearSystem = 
      Teuchos::rcp(new NOX::SharedObject<NOX::Epetra::LinearSystem, 
		                         LOCA::Epetra::Group>(complexLinSys));
  }

  if (finalStatus == NOX::Abstract::Group::Ok)
    isValidComplex = true;

  return finalStatus;

#else

  globalData->locaErrorCheck->throwError(callingFunction, 
					 "Must have EpetraExt support for Hopf tracking.  Configure trilinos with --enable-epetraext");
  return NOX::Abstract::Group::BadDependency;

#endif
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyComplex(const NOX::Abstract::Vector& input_real,
				  const NOX::Abstract::Vector& input_imag,
				  NOX::Abstract::Vector& result_real,
				  NOX::Abstract::Vector& result_imag) const
{
  std::string callingFunction = 
    "LOCA::Epetra::Group::applyComplex()";

  // We must have the time-dependent interface
  if (userInterfaceTime == Teuchos::null)
    return NOX::Abstract::Group::BadDependency;

#ifdef HAVE_NOX_EPETRAEXT
  const NOX::Epetra::Vector& epetra_input_real =
    dynamic_cast<const NOX::Epetra::Vector&>(input_real);
  const NOX::Epetra::Vector& epetra_input_imag =
    dynamic_cast<const NOX::Epetra::Vector&>(input_imag);
  NOX::Epetra::Vector& epetra_result_real =
    dynamic_cast<NOX::Epetra::Vector&>(result_real);
  NOX::Epetra::Vector& epetra_result_imag =
    dynamic_cast<NOX::Epetra::Vector&>(result_imag);

  EpetraExt::BlockVector complex_input(*complexVec);
  complex_input.LoadBlockValues(epetra_input_real.getEpetraVector(), 0);
  complex_input.LoadBlockValues(epetra_input_imag.getEpetraVector(), 1);

  EpetraExt::BlockVector complex_result(*complexVec);
  complex_result.PutScalar(0.0);

  complexMatrix->Apply(complex_input, complex_result);

  complex_result.ExtractBlockValues(epetra_result_real.getEpetraVector(), 0);
  complex_result.ExtractBlockValues(epetra_result_imag.getEpetraVector(), 1);

  return NOX::Abstract::Group::Ok;
#else

  globalData->locaErrorCheck->throwError(callingFunction, 
					 "Must have EpetraExt support for Hopf tracking.  Configure trilinos with --enable-epetraext");
  return NOX::Abstract::Group::BadDependency;

#endif
}

// NOX::Abstract::Group::ReturnType
// LOCA::Epetra::Group::applyComplexMultiVector(
// 				const NOX::Abstract::MultiVector& input_real,
// 				const NOX::Abstract::MultiVector& input_imag,
// 				NOX::Abstract::MultiVector& result_real,
// 				NOX::Abstract::MultiVector& result_imag) const
// {
//   std::string callingFunction = 
//     "LOCA::Epetra::Group::applyComplexMultiVector()";
//   NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
//   NOX::Abstract::Group::ReturnType status;

//   for (int i=0; i<input_real.numVectors(); i++) {
//     status = applyComplex(input_real[i], input_imag[i], 
// 			  result_real[i], result_imag[i]);
//     finalStatus = 
//       globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
// 							     finalStatus,
// 							     callingFunction);
//   }
  
//   return finalStatus;
// }

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyComplexMultiVector(
				const NOX::Abstract::MultiVector& input_real,
				const NOX::Abstract::MultiVector& input_imag,
				NOX::Abstract::MultiVector& result_real,
				NOX::Abstract::MultiVector& result_imag) const
{
  std::string callingFunction = 
    "LOCA::Epetra::Group::applyComplexMultiVector()";
  // We must have the time-dependent interface
  if (userInterfaceTime == Teuchos::null)
    return NOX::Abstract::Group::BadDependency;

#ifdef HAVE_NOX_EPETRAEXT
  const NOX::Epetra::MultiVector& epetra_input_real =
    dynamic_cast<const NOX::Epetra::MultiVector&>(input_real);
  const NOX::Epetra::MultiVector& epetra_input_imag =
    dynamic_cast<const NOX::Epetra::MultiVector&>(input_imag);
  NOX::Epetra::MultiVector& epetra_result_real =
    dynamic_cast<NOX::Epetra::MultiVector&>(result_real);
  NOX::Epetra::MultiVector& epetra_result_imag =
    dynamic_cast<NOX::Epetra::MultiVector&>(result_imag);

  // Get Jacobian matrix for row map
  Teuchos::RCP<Epetra_RowMatrix> jac = Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(sharedLinearSystem.getObject(this)->getJacobianOperator());

  EpetraExt::BlockMultiVector complex_input(jac->RowMatrixRowMap(),
					    complexMatrix->RowMap(),
					    input_real.numVectors());
  complex_input.LoadBlockValues(epetra_input_real.getEpetraMultiVector(), 0);
  complex_input.LoadBlockValues(epetra_input_imag.getEpetraMultiVector(), 1);

  EpetraExt::BlockMultiVector complex_result(jac->RowMatrixRowMap(),
					     complexMatrix->RowMap(),
					     input_real.numVectors());
  complex_result.PutScalar(0.0);

  complexMatrix->Apply(complex_input, complex_result);

  complex_result.ExtractBlockValues(epetra_result_real.getEpetraMultiVector(), 
				    0);
  complex_result.ExtractBlockValues(epetra_result_imag.getEpetraMultiVector(), 
				    1);

  return NOX::Abstract::Group::Ok;
#else

  globalData->locaErrorCheck->throwError(callingFunction, 
					 "Must have EpetraExt support for Hopf tracking.  Configure trilinos with --enable-epetraext");
  return NOX::Abstract::Group::BadDependency;

#endif
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyComplexInverseMultiVector(
			    Teuchos::ParameterList& lsParams,
			    const NOX::Abstract::MultiVector& input_real,
			    const NOX::Abstract::MultiVector& input_imag,
			    NOX::Abstract::MultiVector& result_real,
			    NOX::Abstract::MultiVector& result_imag) const
{
  std::string callingFunction = 
    "LOCA::Epetra::Group::applyComplexInverseMultiVector()";

  // We must have the time-dependent interface
  if (userInterfaceTime == Teuchos::null)
    return NOX::Abstract::Group::BadDependency;

#ifdef HAVE_NOX_EPETRAEXT
  
  // Compute the preconditioner
  NOX::Epetra::LinearSystem::PreconditionerReusePolicyType precPolicy = 
    complexSharedLinearSystem->getObject(this)->getPreconditionerPolicy();

  if (!isValidComplexPrec) {
    if (precPolicy == NOX::Epetra::LinearSystem::PRPT_REBUILD) {
      complexSharedLinearSystem->getObject(this)->destroyPreconditioner();
      complexSharedLinearSystem->getObject(this)->
	createPreconditioner(xVector, 
			     lsParams, 
			     false);
      isValidComplexPrec = true;
    }
    else if (precPolicy == NOX::Epetra::LinearSystem::PRPT_RECOMPUTE) {
      complexSharedLinearSystem->getObject(this)->
	recomputePreconditioner(xVector, lsParams);
    }
    else if (precPolicy == NOX::Epetra::LinearSystem::PRPT_REUSE) {
      // Do Nothing!!!
    }
  }

  EpetraExt::BlockVector complex_input(*complexVec);
  EpetraExt::BlockVector complex_result(*complexVec);
  NOX::Epetra::Vector nev_input(Teuchos::rcp(&complex_input,false),
				NOX::Epetra::Vector::CreateView);
  NOX::Epetra::Vector nev_result(Teuchos::rcp(&complex_result,false), 
				 NOX::Epetra::Vector::CreateView);
  const NOX::Epetra::Vector* epetra_input_real;
  const NOX::Epetra::Vector* epetra_input_imag;
  NOX::Epetra::Vector* epetra_result_real; 
  NOX::Epetra::Vector* epetra_result_imag; 
  bool status;
  bool finalStatus = true;
  for (int i=0; i<input_real.numVectors(); i++) {
    epetra_input_real = 
      dynamic_cast<const NOX::Epetra::Vector*>(&input_real[i]);
    epetra_input_imag = 
      dynamic_cast<const NOX::Epetra::Vector*>(&input_imag[i]);
    epetra_result_real = dynamic_cast<NOX::Epetra::Vector*>(&result_real[i]);
    epetra_result_imag = dynamic_cast<NOX::Epetra::Vector*>(&result_imag[i]);

    complex_input.LoadBlockValues(epetra_input_real->getEpetraVector(), 0);
    complex_input.LoadBlockValues(epetra_input_imag->getEpetraVector(), 1);
    complex_result.PutScalar(0.0);
    
    status = 
      complexSharedLinearSystem->getObject(this)->applyJacobianInverse(
							 lsParams, 
							 nev_input,
							 nev_result);

    complex_result.ExtractBlockValues(epetra_result_real->getEpetraVector(),0);
    complex_result.ExtractBlockValues(epetra_result_imag->getEpetraVector(),1);
    
    finalStatus = finalStatus && status;
  }
  
  if (finalStatus)
    return NOX::Abstract::Group::Ok;
  else
    return NOX::Abstract::Group::NotConverged;
#else

  globalData->locaErrorCheck->throwError(callingFunction, 
					 "Must have EpetraExt support for Hopf tracking.  Configure trilinos with --enable-epetraext");
  return NOX::Abstract::Group::BadDependency;

#endif
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyComplexTranspose(
				  const NOX::Abstract::Vector& input_real,
				  const NOX::Abstract::Vector& input_imag,
				  NOX::Abstract::Vector& result_real,
				  NOX::Abstract::Vector& result_imag) const
{
  std::string callingFunction = 
    "LOCA::Epetra::Group::applyComplexTranspose()";

  // We must have the time-dependent interface
  if (userInterfaceTime == Teuchos::null)
    return NOX::Abstract::Group::BadDependency;

#ifdef HAVE_NOX_EPETRAEXT
  const NOX::Epetra::Vector& epetra_input_real =
    dynamic_cast<const NOX::Epetra::Vector&>(input_real);
  const NOX::Epetra::Vector& epetra_input_imag =
    dynamic_cast<const NOX::Epetra::Vector&>(input_imag);
  NOX::Epetra::Vector& epetra_result_real =
    dynamic_cast<NOX::Epetra::Vector&>(result_real);
  NOX::Epetra::Vector& epetra_result_imag =
    dynamic_cast<NOX::Epetra::Vector&>(result_imag);

  EpetraExt::BlockVector complex_input(*complexVec);
  complex_input.LoadBlockValues(epetra_input_real.getEpetraVector(), 0);
  complex_input.LoadBlockValues(epetra_input_imag.getEpetraVector(), 1);

  EpetraExt::BlockVector complex_result(*complexVec);
  complex_result.PutScalar(0.0);

  // Note:  transpose of block matrix is equivalent to conjugate transpose
  // of complex matrix.  Hence we don't need to take a conjugate

  complexMatrix->SetUseTranspose(true);
  complexMatrix->Apply(complex_input, complex_result);
  complexMatrix->SetUseTranspose(false);

  complex_result.ExtractBlockValues(epetra_result_real.getEpetraVector(), 0);
  complex_result.ExtractBlockValues(epetra_result_imag.getEpetraVector(), 1);
  
  return NOX::Abstract::Group::Ok;
#else

  globalData->locaErrorCheck->throwError(callingFunction, 
					 "Must have EpetraExt support for Hopf tracking.  Configure trilinos with --enable-epetraext");
  return NOX::Abstract::Group::BadDependency;

#endif
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyComplexTransposeMultiVector(
				const NOX::Abstract::MultiVector& input_real,
				const NOX::Abstract::MultiVector& input_imag,
				NOX::Abstract::MultiVector& result_real,
				NOX::Abstract::MultiVector& result_imag) const
{
  std::string callingFunction = 
    "LOCA::Epetra::Group::applyComplexTransposeMultiVector()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  for (int i=0; i<input_real.numVectors(); i++) {
    status = applyComplexTranspose(input_real[i], input_imag[i], 
				   result_real[i], result_imag[i]);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }
  
  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyComplexTransposeInverseMultiVector(
			    Teuchos::ParameterList& lsParams,
			    const NOX::Abstract::MultiVector& input_real,
			    const NOX::Abstract::MultiVector& input_imag,
			    NOX::Abstract::MultiVector& result_real,
			    NOX::Abstract::MultiVector& result_imag) const
{
  std::string callingFunction = 
    "LOCA::Epetra::Group::applyComplexTransposeInverseMultiVector()";

  // We must have the time-dependent interface
  if (userInterfaceTime == Teuchos::null)
    return NOX::Abstract::Group::BadDependency;

#ifdef HAVE_NOX_EPETRAEXT
  if (complexSharedLinearSystem == Teuchos::null) {

    NOX::Epetra::Vector nev(complexVec, NOX::Epetra::Vector::CreateView);

    // Create the Linear System
    Teuchos::RCP<NOX::Epetra::Interface::Required> iReq;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac;
    Teuchos::RCP<NOX::Epetra::LinearSystem> complexLinSys = 
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, 
							lsParams,
							iReq, 
							iJac, 
							complexMatrix, 
							nev));
    complexLinSys->setJacobianOperatorForSolve(complexMatrix);
    complexSharedLinearSystem = 
      Teuchos::rcp(new NOX::SharedObject<NOX::Epetra::LinearSystem, 
		                         LOCA::Epetra::Group>(complexLinSys));
  }
  Teuchos::RCP<NOX::Epetra::LinearSystem> complexLinSys = 
    Teuchos::rcp_const_cast<NOX::Epetra::LinearSystem>(complexSharedLinearSystem->getObject());

  // Instantiate transpose solver
  LOCA::Epetra::TransposeLinearSystem::Factory tls_factory(globalData);
  if (tls_strategy == Teuchos::null)
    tls_strategy = tls_factory.create(Teuchos::rcp(&lsParams, false), 
				      complexLinSys);

  // Compute complex transpose
  tls_strategy->createJacobianTranspose();

  // Now compute preconditioner
  tls_strategy->createTransposePreconditioner(xVector, lsParams);

  // Solve each complex system
  EpetraExt::BlockVector complex_input(*complexVec);
  EpetraExt::BlockVector complex_result(*complexVec);
  NOX::Epetra::Vector nev_input(Teuchos::rcp(&complex_input,false),
				NOX::Epetra::Vector::CreateView);
  NOX::Epetra::Vector nev_result(Teuchos::rcp(&complex_result,false), 
				 NOX::Epetra::Vector::CreateView);
  const NOX::Epetra::Vector* epetra_input_real;
  const NOX::Epetra::Vector* epetra_input_imag;
  NOX::Epetra::Vector* epetra_result_real; 
  NOX::Epetra::Vector* epetra_result_imag; 
  bool status;
  bool finalStatus = true;
  for (int i=0; i<input_real.numVectors(); i++) {
    epetra_input_real = 
      dynamic_cast<const NOX::Epetra::Vector*>(&input_real[i]);
    epetra_input_imag = 
      dynamic_cast<const NOX::Epetra::Vector*>(&input_imag[i]);
    epetra_result_real = dynamic_cast<NOX::Epetra::Vector*>(&result_real[i]);
    epetra_result_imag = dynamic_cast<NOX::Epetra::Vector*>(&result_imag[i]);

    complex_input.LoadBlockValues(epetra_input_real->getEpetraVector(), 0);
    complex_input.LoadBlockValues(epetra_input_imag->getEpetraVector(), 1);
    complex_result.PutScalar(0.0);
    
    status = 
      tls_strategy->applyJacobianTransposeInverse(lsParams, 
						  nev_input,
						  nev_result);

    complex_result.ExtractBlockValues(epetra_result_real->getEpetraVector(),0);
    complex_result.ExtractBlockValues(epetra_result_imag->getEpetraVector(),1);
    
    finalStatus = finalStatus && status;
  }

  // Set original operators in linear system
  complexMatrix->SetUseTranspose(false);
  complexLinSys->setJacobianOperatorForSolve(complexMatrix);
  complexLinSys->destroyPreconditioner();
  
  if (finalStatus)
    return NOX::Abstract::Group::Ok;
  else
    return NOX::Abstract::Group::NotConverged;
#else

  globalData->locaErrorCheck->throwError(callingFunction, 
					 "Must have EpetraExt support for Hopf tracking.  Configure trilinos with --enable-epetraext");
  return NOX::Abstract::Group::BadDependency;

#endif
}

// NOX::Abstract::Group::ReturnType
// LOCA::Epetra::Group::applyShiftedMatrixInverse(
// 					    Teuchos::ParameterList& params,
// 					    const NOX::Abstract::Vector& input,
// 					    NOX::Abstract::Vector& result,
// 					    double shift)
                                           
// {
//   const NOX::Epetra::Vector& epetraInput = 
//     dynamic_cast<const NOX::Epetra::Vector&>(input);
//   NOX::Epetra::Vector& epetraResult = 
//     dynamic_cast<NOX::Epetra::Vector&>(result);

//   // If shift is zero, just apply Jacobian inverse

//   if(shift == 0.0) {
//     applyJacobianInverse(params, epetraInput, epetraResult);
//   }
//   else {

//     // Otherwise, construct a shift and invert operator, and use AztecOO to 
//     // solve linear system

//     Teuchos::RCP<LOCA::Epetra::ShiftInvertOperator> A =
//       Teuchos::rcp(new LOCA::Epetra::ShiftInvertOperator(
// 		  globalData,
// 		  Teuchos::rcp(this,false),
// 		  sharedLinearSystem.getObject(this)->getJacobianOperator(),
// 		  shift,
// 		  userInterfaceTime != Teuchos::null));

//     NOX::Epetra::Vector dummy(epetraResult, NOX::ShapeCopy);
//     Epetra_Vector& epetra_dummy = dummy.getEpetraVector();    
//     Teuchos::RCP<LOCA::Epetra::ShiftInvertInterface> interface = 
//       Teuchos::rcp(new LOCA::Epetra::ShiftInvertInterface); 
//     Teuchos::ParameterList& solveList = params.sublist("NOX").sublist("Direction").sublist("Newton").sublist("Linear Solver");

//     NOX::Epetra::LinearSystemAztecOO shiftsys(
// 	params,
// 	solveList,
// 	userInterface,
// 	Teuchos::rcp_dynamic_cast<NOX::Epetra::Interface::Jacobian>(interface),
// 	A,
// 	epetra_dummy); 

//     shiftsys.setJacobianOperatorForSolve(A);

//     shiftsys.applyJacobianInverse(params,epetraInput,epetraResult);

//   }

//   return NOX::Abstract::Group::Ok;  
// }

Teuchos::RCP<NOX::Epetra::Interface::Required>
LOCA::Epetra::Group::getUserInterface()
{
  return userInterface;
}

void
LOCA::Epetra::Group::printSolution(const NOX::Epetra::Vector& x_,
				      const double conParam) const
{
  userInterface->printSolution(x_.getEpetraVector(), conParam);
}

void
LOCA::Epetra::Group::setScaleVector(const NOX::Abstract::Vector& s)
{
  scaleVecPtr = s.clone(NOX::DeepCopy);

  return;
}

void
LOCA::Epetra::Group::setJacobianOperatorForSolve(
		  const Teuchos::RCP<const Epetra_Operator>& op) const
{
  // Set Jacobian operator for solve
  sharedLinearSystem.getObject(this)->setJacobianOperatorForSolve(op);
  isValidSolverJacOp = true;
}

void
LOCA::Epetra::Group::setFreeEnergyInterface(
              const Teuchos::RCP<LOCA::Epetra::Interface::FreeEnergy>& iFE)
{
  userInterfaceFreeEnergy = iFE;
}

void
LOCA::Epetra::Group::declareSeparateMatrixMemory(bool separateMem)
{
  separateMatrixMemoryDeclared = separateMem;

  // Error check
  if ( separateMem) {
    Teuchos::RCP<Epetra_Operator> linsysjac = 
      sharedLinearSystem.getObject(this)->getJacobianOperator();
    Teuchos::RCP<Epetra_Operator> shiftedlinsysjac = 
      shiftedSharedLinearSystem->getObject(this)->getJacobianOperator();
  
    if (shiftedlinsysjac.get() == linsysjac.get())
      globalData->locaErrorCheck->throwError("LOCA::Epetra::Group::declareSeparateMatrixMemory", 
  	   "User code has called declareSeparateMatrixMemory(), but the pointers to "
           "the matrices in both\nlinear systems are the same! Need to allocate separate memory "
           "in linSys and shiftedlinSys.");
  }
}

double
LOCA::Epetra::Group::computeFreeEnergy()
{
  // We store a real shifted matrix
  if (userInterfaceFreeEnergy != Teuchos::null) {

     userInterface->setParameters(params);

     return userInterfaceFreeEnergy->computeFreeEnergy(xVector.getEpetraVector());
  }
  else {
    globalData->locaErrorCheck->throwError("LOCA::Epetra::Group::computeFreeEnergy", 
	 "Free Energy Calculation needed by Phase Transition algorithm."
         "LOCA::Epetra::Group needs LOCA::Epetra::Interface::FreeEnergy");
  }
  return 0; // not reached
}

void
LOCA::Epetra::Group::resetIsValid()
{
  NOX::Epetra::Group::resetIsValid();
  isValidShiftedPrec = false;
  isValidComplex = false;
  isValidComplexPrec = false;
}

Teuchos::RCP<const NOX::Epetra::LinearSystem>
LOCA::Epetra::Group::getComplexLinearSystem() const
{
  return complexSharedLinearSystem->getObject();
}

Teuchos::RCP<NOX::Epetra::LinearSystem>
LOCA::Epetra::Group::getComplexLinearSystem()
{
  return complexSharedLinearSystem->getObject(this);
}

void
LOCA::Epetra::Group::getComplexMaps(
		 Teuchos::RCP<const Epetra_BlockMap>& baseMap,
		 Teuchos::RCP<const Epetra_BlockMap>& globalMap) const
{
#ifdef HAVE_NOX_EPETRAEXT
  baseMap = Teuchos::rcp(&(xVector.getEpetraVector().Map()),false);
  globalMap = Teuchos::rcp(&(complexVec->Map()),false);
#else
  globalData->locaErrorCheck->throwError("LOCA::Epetra::Group::getComplexMaps()", 
					 "Must have EpetraExt support for Hopf tracking.  Configure trilinos with --enable-epetraext");
#endif
}
