// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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

#ifdef HAVE_NOX_EPETRAEXT
#include "EpetraExt_BlockCrsMatrix.h"
#include "EpetraExt_BlockVector.h"
#endif

LOCA::Epetra::Group::Group(
	    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	    Teuchos::ParameterList& printingParams, 
	    const Teuchos::RefCountPtr<LOCA::Epetra::Interface::Required>& i, 
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
  isValidComplexPrec(false)
{
}

LOCA::Epetra::Group::Group(
	    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	    Teuchos::ParameterList& printingParams, 
	    const Teuchos::RefCountPtr<LOCA::Epetra::Interface::Required>& i, 
	    NOX::Epetra::Vector& initialGuess, 
	    const Teuchos::RefCountPtr<NOX::Epetra::LinearSystem>& linSys,
	    const LOCA::ParameterVector& p) :
  NOX::Epetra::Group(printingParams, i, initialGuess, linSys),
  LOCA::Abstract::Group(global_data),
  globalData(global_data),
  printParams(printingParams),
  params(p),
  userInterface(i),
  userInterfaceTime(),
  userInterfaceTimeMF(),
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
  isValidComplexPrec(false)
{
}

LOCA::Epetra::Group::Group(
	const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	Teuchos::ParameterList& printingParams, 
	const Teuchos::RefCountPtr<LOCA::Epetra::Interface::TimeDependent>& i, 
	NOX::Epetra::Vector& initialGuess, 
	const Teuchos::RefCountPtr<NOX::Epetra::LinearSystem>& linSys,
	const Teuchos::RefCountPtr<NOX::Epetra::LinearSystem>& shiftedLinSys,
	const LOCA::ParameterVector& p) :
  NOX::Epetra::Group(printingParams, i, initialGuess, linSys),
  LOCA::Abstract::Group(global_data),
  globalData(global_data),
  printParams(printingParams),
  params(p),
  userInterface(i), 
  userInterfaceTime(i),
  userInterfaceTimeMF(),
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
  isValidComplexPrec(false)
{
  shiftedSharedLinearSystem = 
    Teuchos::rcp(new NOX::SharedObject<NOX::Epetra::LinearSystem, 
		                       LOCA::Epetra::Group>(shiftedLinSys));
}

LOCA::Epetra::Group::Group(
	    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	    Teuchos::ParameterList& printingParams, 
	    const Teuchos::RefCountPtr<LOCA::Epetra::Interface::TimeDependentMatrixFree>& i, 
	    NOX::Epetra::Vector& initialGuess, 
	    const Teuchos::RefCountPtr<NOX::Epetra::LinearSystem>& linSys,
	    const Teuchos::RefCountPtr<NOX::Epetra::LinearSystem>& shiftedLinSys,
	    const LOCA::ParameterVector& p) :
  NOX::Epetra::Group(printingParams, i, initialGuess, linSys),
  LOCA::Abstract::Group(global_data),
  globalData(global_data),
  printParams(printingParams),
  params(p),
  userInterface(i),
  userInterfaceTime(),
  userInterfaceTimeMF(i),
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
  isValidComplexPrec(false)
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
  isValidComplexPrec(source.isValidComplexPrec)
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

Teuchos::RefCountPtr<NOX::Abstract::Group>
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
  string callingFunction = 
    "LOCA::Epetra::Group::applyJacobianTransposeInverse()";
  NOX::Abstract::Group::ReturnType status;
  
  // Get non-const linsys
  Teuchos::RefCountPtr<NOX::Epetra::LinearSystem> linSys = 
    sharedLinearSystem.getObject(this);

  // Get Jacobian operator
  Teuchos::RefCountPtr<Epetra_Operator> jac =
    linSys->getJacobianOperator();

  // Instantiate transpose solver
  LOCA::Epetra::TransposeLinearSystem::Factory tls_factory(globalData);
  if (tls_strategy == Teuchos::null)
    tls_strategy = tls_factory.create(Teuchos::rcp(&params, false), linSys);
     
  // Compute Jacobian transpose J^T
  tls_strategy->computeJacobianTranspose(xVector);

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
  string callingFunction = 
    "LOCA::Epetra::Group::applyJacobianTransposeInverseMultiVector()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // Get non-const linsys
  Teuchos::RefCountPtr<NOX::Epetra::LinearSystem> linSys = 
    sharedLinearSystem.getObject(this);

  // Get Jacobian operator
  Teuchos::RefCountPtr<Epetra_Operator> jac =
    linSys->getJacobianOperator();

  // Instantiate transpose solver
  LOCA::Epetra::TransposeLinearSystem::Factory tls_factory(globalData);
  if (tls_strategy == Teuchos::null)
    tls_strategy = tls_factory.create(Teuchos::rcp(&params, false), linSys);
     
  // Compute Jacobian transpose J^T
  tls_strategy->computeJacobianTranspose(xVector);

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
LOCA::Epetra::Group::setParam(string paramID, double val)
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
LOCA::Epetra::Group::getParam(string paramID) const
{
  return params.getValue(paramID);
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
    Teuchos::RefCountPtr<NOX::Abstract::Vector> as = a.clone(NOX::DeepCopy);
    Teuchos::RefCountPtr<NOX::Abstract::Vector> bs = b.clone(NOX::DeepCopy);
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
LOCA::Epetra::Group::augmentJacobianForHomotopy(double p)
{

  //Allocate temporary vectors if not aready done
  if (tmpVectorPtr == Teuchos::null)
    tmpVectorPtr = Teuchos::rcp(new Epetra_Vector(xVector.getEpetraVector()));
  if (tmpVectorPtr2 == Teuchos::null)
    tmpVectorPtr2 = Teuchos::rcp(new Epetra_Vector(xVector.getEpetraVector()));

  tmpVectorPtr2->PutScalar(1.0-p);

  // See if it is an Epetra_CrsMatrix
  Teuchos::RefCountPtr<const Epetra_CrsMatrix> constTestCrs;
  Teuchos::RefCountPtr<Epetra_CrsMatrix> testCrs;
  constTestCrs = Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>
    (sharedLinearSystem.getObject(this)->getJacobianOperator());
  if (constTestCrs != Teuchos::null) {
    testCrs = Teuchos::rcp_const_cast<Epetra_CrsMatrix>(constTestCrs);
    testCrs->Scale(p);
    testCrs->ExtractDiagonalCopy(*tmpVectorPtr);
    tmpVectorPtr->Update(1.0, *tmpVectorPtr2, 1.0);
    testCrs->ReplaceDiagonalValues(*tmpVectorPtr);
    return LOCA::Abstract::Group::Ok;

  }

  // See if it is an Epetra_VbrMatrix
  Teuchos::RefCountPtr<const Epetra_VbrMatrix> constTestVbr;
  Teuchos::RefCountPtr<Epetra_VbrMatrix> testVbr;
  constTestVbr = Teuchos::rcp_dynamic_cast<const Epetra_VbrMatrix>
    (sharedLinearSystem.getObject(this)->getJacobianOperator());
  if (constTestVbr != Teuchos::null) {
    testVbr = Teuchos::rcp_const_cast<Epetra_VbrMatrix>(constTestVbr);
    testVbr->Scale(p);
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
    Teuchos::RefCountPtr<Epetra_Operator> mass = 
      shiftedSharedLinearSystem->getObject(this)->getJacobianOperator();

    bool res = 
      userInterfaceTime->computeShiftedMatrix(alpha, beta, 
					      xVector.getEpetraVector(), 
					      *mass);
    
    // Check if Jacobian and mass matrices are the same, in which case
    // the Jacobian is no longer valid
    Teuchos::RefCountPtr<Epetra_Operator> jac = 
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
  string callingFunction = 
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
    bool reusePrec = 
      shiftedSharedLinearSystem->getObject(this)->checkPreconditionerReuse();

    if (!isValidShiftedPrec  && !reusePrec ) {
      shiftedSharedLinearSystem->getObject(this)->destroyPreconditioner();
      shiftedSharedLinearSystem->getObject(this)->
	createPreconditioner(xVector, lsParams, false);
      isValidShiftedPrec = true;
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

bool
LOCA::Epetra::Group::isComplex() const
{
  return isValidComplex;
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::computeComplex(double frequency)
{
  string callingFunction = 
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
  Teuchos::RefCountPtr<Epetra_RowMatrix> jac = Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(sharedLinearSystem.getObject(this)->getJacobianOperator());
   
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
  Teuchos::RefCountPtr<Epetra_RowMatrix> mass = Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(shiftedSharedLinearSystem->getObject(this)->getJacobianOperator());

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
  string callingFunction = 
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

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyComplexMultiVector(
				const NOX::Abstract::MultiVector& input_real,
				const NOX::Abstract::MultiVector& input_imag,
				NOX::Abstract::MultiVector& result_real,
				NOX::Abstract::MultiVector& result_imag) const
{
  string callingFunction = 
    "LOCA::Epetra::Group::applyComplexMultiVector()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  for (int i=0; i<input_real.numVectors(); i++) {
    status = applyComplex(input_real[i], input_imag[i], 
			  result_real[i], result_imag[i]);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }
  
  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyComplexInverseMultiVector(
			    Teuchos::ParameterList& lsParams,
			    const NOX::Abstract::MultiVector& input_real,
			    const NOX::Abstract::MultiVector& input_imag,
			    NOX::Abstract::MultiVector& result_real,
			    NOX::Abstract::MultiVector& result_imag) const
{
  string callingFunction = 
    "LOCA::Epetra::Group::applyComplexInverseMultiVector()";

  // We must have the time-dependent interface
  if (userInterfaceTime == Teuchos::null)
    return NOX::Abstract::Group::BadDependency;

#ifdef HAVE_NOX_EPETRAEXT
  if (complexSharedLinearSystem == Teuchos::null) {

    NOX::Epetra::Vector nev(complexVec, NOX::Epetra::Vector::CreateView);

    // Create the Linear System
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Required> iReq;
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> iJac;
    Teuchos::RefCountPtr<NOX::Epetra::LinearSystem> complexLinSys = 
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

  bool reusePrec = 
      complexSharedLinearSystem->getObject(this)->checkPreconditionerReuse();

  if (!isValidComplexPrec  && !reusePrec ) {
    complexSharedLinearSystem->getObject(this)->destroyPreconditioner();
    complexSharedLinearSystem->getObject(this)->createPreconditioner(xVector, 
								    lsParams, 
								    false);
    isValidComplexPrec = true;
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

//     Teuchos::RefCountPtr<LOCA::Epetra::ShiftInvertOperator> A =
//       Teuchos::rcp(new LOCA::Epetra::ShiftInvertOperator(
// 		  globalData,
// 		  Teuchos::rcp(this,false),
// 		  sharedLinearSystem.getObject(this)->getJacobianOperator(),
// 		  shift,
// 		  userInterfaceTime != Teuchos::null));

//     NOX::Epetra::Vector dummy(epetraResult, NOX::ShapeCopy);
//     Epetra_Vector& epetra_dummy = dummy.getEpetraVector();    
//     Teuchos::RefCountPtr<LOCA::Epetra::ShiftInvertInterface> interface = 
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

Teuchos::RefCountPtr<NOX::Epetra::Interface::Required>
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
		  const Teuchos::RefCountPtr<const Epetra_Operator>& op) const
{
  // Set Jacobian operator for solve
  sharedLinearSystem.getObject(this)->setJacobianOperatorForSolve(op);
  isValidSolverJacOp = true;
}

void
LOCA::Epetra::Group::resetIsValid()
{
  NOX::Epetra::Group::resetIsValid();
  isValidShiftedPrec = false;
  isValidComplex = false;
  isValidComplexPrec = false;
}
