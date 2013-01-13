// $Id$
// $Source$

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

#include "LOCA_Hopf_MooreSpence_ExtendedGroup.H"
#include "LOCA_Hopf_MooreSpence_AbstractGroup.H"
#include "LOCA_Hopf_MooreSpence_SolverStrategy.H"
#include "LOCA_Parameter_Vector.H"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_GlobalData.H"
#include "LOCA_Factory.H"
#include "LOCA_Parameter_SublistParser.H"
#include "NOX_Utils.H"
#include "LOCA_ErrorCheck.H"

LOCA::Hopf::MooreSpence::ExtendedGroup::ExtendedGroup(
	 const Teuchos::RCP<LOCA::GlobalData>& global_data,
         const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& hpfParams,
	 const Teuchos::RCP<LOCA::Hopf::MooreSpence::AbstractGroup>& g)
  : LOCA::Extended::MultiAbstractGroup(),
    LOCA::MultiContinuation::AbstractGroup(),
    globalData(global_data),
    parsedParams(topParams),
    hopfParams(hpfParams),
    grpPtr(g),
    xMultiVec(globalData, g->getX(), 1),
    fMultiVec(globalData, g->getX(), 2),
    newtonMultiVec(globalData, g->getX(), 1),
    lengthMultiVec(),
    xVec(),
    fVec(),
    ffMultiVec(),
    dfdpMultiVec(),
    newtonVec(),
    lengthVec(),
    massTimesY(),
    minusMassTimesZ(),
    solverStrategy(),
    index_f(1),
    index_dfdp(1),
    bifParamID(1), 
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false)
{
  const char *func = "LOCA::Hopf::MooreSpence::ExtendedGroup()";

  // Set x
  *(xMultiVec.getColumn(0)->getXVec()) = g->getX();

  if (!hopfParams->isParameter("Bifurcation Parameter")) {
    globalData->locaErrorCheck->throwError(func,
				 "\"Bifurcation Parameter\" name is not set!");
  }
  std::string bifParamName = hopfParams->get("Bifurcation Parameter", "None");
  const ParameterVector& p = grpPtr->getParams();
  bifParamID[0] = p.getIndex(bifParamName);

  if (!hopfParams->isParameter("Length Normalization Vector")) {
    globalData->locaErrorCheck->throwError(func,
			   "\"Length Normalization Vector\" is not set!");
  }
  Teuchos::RCP<NOX::Abstract::Vector> lenVecPtr = 
    (*hopfParams).INVALID_TEMPLATE_QUALIFIER 
    get< Teuchos::RCP<NOX::Abstract::Vector> >(
					       "Length Normalization Vector");

  if (!hopfParams->isParameter("Initial Real Eigenvector")) {
    globalData->locaErrorCheck->throwError(func,
				 "\"Initial Real Eigenvector\" is not set!");
  }
  Teuchos::RCP<NOX::Abstract::Vector> realEigVecPtr = 
    (*hopfParams).INVALID_TEMPLATE_QUALIFIER 
    get<Teuchos::RCP<NOX::Abstract::Vector> >(
						  "Initial Real Eigenvector");

  if (!hopfParams->isParameter("Initial Imaginary Eigenvector")) {
    globalData->locaErrorCheck->throwError(func,
			      "\"Initial Imaginary Eigenvector\" is not set!");
  }
  Teuchos::RCP<NOX::Abstract::Vector> imagEigVecPtr = 
    (*hopfParams).INVALID_TEMPLATE_QUALIFIER 
    get<Teuchos::RCP<NOX::Abstract::Vector> >(
					       "Initial Imaginary Eigenvector");

  if (!hopfParams->isParameter("Initial Frequency")) {
    globalData->locaErrorCheck->throwError(func,
				 "\"Initial Frequency\" is not set!");
  }
  xMultiVec.getColumn(0)->getFrequency() = 
    hopfParams->get("Initial Frequency", 1.0);

  bool perturbSoln = hopfParams->get("Perturb Initial Solution", false);
  double perturbSize = hopfParams->get("Relative Perturbation Size", 1.0e-3);

  lengthMultiVec = 
    lenVecPtr->createMultiVector(1, NOX::DeepCopy);
  *(xMultiVec.getColumn(0)->getRealEigenVec()) = *realEigVecPtr;
  *(xMultiVec.getColumn(0)->getImagEigenVec()) = *imagEigVecPtr;
  massTimesY = lengthMultiVec->clone(NOX::ShapeCopy);
  minusMassTimesZ = lengthMultiVec->clone(NOX::ShapeCopy);

  // Instantiate solver strategy
  solverStrategy = 
    globalData->locaFactory->createMooreSpenceHopfSolverStrategy(parsedParams,
								 hopfParams);

  // Set up multi-vector views
  setupViews(); 

  init(perturbSoln, perturbSize);
}

LOCA::Hopf::MooreSpence::ExtendedGroup::ExtendedGroup(
		const LOCA::Hopf::MooreSpence::ExtendedGroup& source, 
		NOX::CopyType type)
  : globalData(source.globalData),
    parsedParams(source.parsedParams),
    hopfParams(source.hopfParams),
    grpPtr(Teuchos::rcp_dynamic_cast<LOCA::Hopf::MooreSpence::AbstractGroup>(source.grpPtr->clone(type))),
    xMultiVec(source.xMultiVec, type),
    fMultiVec(source.fMultiVec, type),
    newtonMultiVec(source.newtonMultiVec, type),
    lengthMultiVec(source.lengthMultiVec->clone(type)),
    xVec(),
    fVec(),
    ffMultiVec(),
    dfdpMultiVec(),
    newtonVec(),
    lengthVec(),
    massTimesY(source.massTimesY->clone(type)),
    minusMassTimesZ(source.minusMassTimesZ->clone(type)),
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
    globalData->locaFactory->createMooreSpenceHopfSolverStrategy(parsedParams,
								 hopfParams);

  // Set up multi-vector views
  setupViews();

  if (type == NOX::ShapeCopy) {
    isValidF = false;
    isValidJacobian = false;
    isValidNewton = false;
  }
}

LOCA::Hopf::MooreSpence::ExtendedGroup::~ExtendedGroup() 
{
}

NOX::Abstract::Group&
LOCA::Hopf::MooreSpence::ExtendedGroup::operator=(
					   const NOX::Abstract::Group& source)
{
  copy(source);
  return *this;
}

Teuchos::RCP<NOX::Abstract::Group>
LOCA::Hopf::MooreSpence::ExtendedGroup::clone(NOX::CopyType type) const 
{
  return 
    Teuchos::rcp(new LOCA::Hopf::MooreSpence::ExtendedGroup(*this, type));
}

void
LOCA::Hopf::MooreSpence::ExtendedGroup::setX(const NOX::Abstract::Vector& y) 
{
  const LOCA::Hopf::MooreSpence::ExtendedVector& yy = 
    dynamic_cast<const LOCA::Hopf::MooreSpence::ExtendedVector&>(y);
  grpPtr->setX( *yy.getXVec() );
  *xVec = y;
  setBifParam(xVec->getBifParam());

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

void
LOCA::Hopf::MooreSpence::ExtendedGroup::computeX(
					      const NOX::Abstract::Group& g, 
					      const NOX::Abstract::Vector& d,
					      double step) 
{
  const LOCA::Hopf::MooreSpence::ExtendedGroup& gg = 
    dynamic_cast<const LOCA::Hopf::MooreSpence::ExtendedGroup&>(g);
  const LOCA::Hopf::MooreSpence::ExtendedVector& dd = 
    dynamic_cast<const LOCA::Hopf::MooreSpence::ExtendedVector&>(d);

  grpPtr->computeX(*(gg.grpPtr), *dd.getXVec(), step);
  xVec->update(1.0, gg.getX(), step, dd, 0.0);
  setBifParam(xVec->getBifParam());

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::computeF() 
{
  if (isValidF)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction = 
    "LOCA::Hopf::MooreSpence::ExtendedGroup::computeF()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Compute underlying F
  if (!grpPtr->isF()) {
    status = grpPtr->computeF();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }
  *(fVec->getXVec()) = grpPtr->getF();
  
  // Compute underlying (J+iwB)
  if (!grpPtr->isComplex()) {
    status = grpPtr->computeComplex(xVec->getFrequency());
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  // Compute (J+iwB)*(y+iz)
  status = grpPtr->applyComplex(*(xVec->getRealEigenVec()), 
				*(xVec->getImagEigenVec()),
				*(fVec->getRealEigenVec()),
				*(fVec->getImagEigenVec()));
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  
  // Compute l^T*y - 1
  fVec->getFrequency() = lTransNorm(*(xVec->getRealEigenVec())) - 1.0;

  // Compute l^T*z - 1
  fVec->getBifParam() = lTransNorm(*(xVec->getImagEigenVec()));
  
  isValidF = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::computeJacobian() 
{
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction = 
    "LOCA::Hopf::MooreSpence::ExtendedGroup::computeJacobian()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Compute underlying df/dp (may invalidate underlying data)
  status = grpPtr->computeDfDpMulti(bifParamID, 
				    *fMultiVec.getXMultiVec(), 
				    isValidF);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute underlying d(J+iwB)*(y+iz)/dp (may invalidate underlying data)
  status = grpPtr->computeDCeDp(bifParamID,
				*(xVec->getRealEigenVec()), 
				*(xVec->getImagEigenVec()),
				xVec->getFrequency(),
				*(fMultiVec.getRealEigenMultiVec()),
				*(fMultiVec.getImagEigenMultiVec()),
				isValidF);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  // Compute underlying Jacobian
  status = grpPtr->computeJacobian();
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  // Compute underlying mass matrix
  status = grpPtr->computeShiftedMatrix(0.0, 1.0);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  // Compute B*y
  status = 
    grpPtr->applyShiftedMatrixMultiVector(*(xMultiVec.getRealEigenMultiVec()),
					  *massTimesY);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  // Compute -B*z
  status = 
    grpPtr->applyShiftedMatrixMultiVector(*(xMultiVec.getImagEigenMultiVec()),
					  *minusMassTimesZ);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  minusMassTimesZ->scale(-1.0);

  // Compute complex matrix
  status = grpPtr->computeComplex(xVec->getFrequency());
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  solverStrategy->setBlocks(grpPtr, 
			    Teuchos::rcp(this, false),
			    xVec->getRealEigenVec(),
			    xVec->getImagEigenVec(),
			    fVec->getRealEigenVec(),
			    fVec->getImagEigenVec(),
			    fMultiVec.getColumn(1)->getXVec(),
			    fMultiVec.getColumn(1)->getRealEigenVec(),
			    fMultiVec.getColumn(1)->getImagEigenVec(),
			    Teuchos::rcp(&((*massTimesY)[0]),false),
			    Teuchos::rcp(&((*minusMassTimesZ)[0]),false),
			    xVec->getFrequency());

  isValidJacobian = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::computeGradient() 
{
  return NOX::Abstract::Group::NotDefined;
}
   
NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::computeNewton(
						 Teuchos::ParameterList& params) 
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction = 
    "LOCA::Hopf::MooreSpence::ExtendedGroup::computeNewton()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Make sure F is valid
  if (!isF()) {
    status = computeF();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }
  
  // Make sure Jacobian is valid
  if (!isJacobian()) {
    status = computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  // zero out newton vec -- used as initial guess for some linear solvers
  newtonMultiVec.init(0.0);

  // solve
  status = solverStrategy->solve(params, *ffMultiVec, newtonMultiVec);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  newtonMultiVec.scale(-1.0);

  isValidNewton = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::applyJacobian(
					  const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  // Convert input, result to multivectors
  Teuchos::RCP<NOX::Abstract::MultiVector> mv_input = 
    input.createMultiVector(1, NOX::DeepCopy);
  Teuchos::RCP<NOX::Abstract::MultiVector> mv_result = 
    result.createMultiVector(1, NOX::DeepCopy);

  // Call multivector version of applyJacobian
  NOX::Abstract::Group::ReturnType status = 
    applyJacobianMultiVector(*mv_input, *mv_result);

  // Copy result
  result = (*mv_result)[0];

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::applyJacobianTranspose(
					  const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  // Convert input, result to multivectors
  Teuchos::RCP<NOX::Abstract::MultiVector> mv_input = 
    input.createMultiVector(1, NOX::DeepCopy);
  Teuchos::RCP<NOX::Abstract::MultiVector> mv_result = 
    result.createMultiVector(1, NOX::DeepCopy);

  // Call multivector version of applyJacobianTranspose
  NOX::Abstract::Group::ReturnType status = 
    applyJacobianTransposeMultiVector(*mv_input, *mv_result);

  // Copy result
  result = (*mv_result)[0];

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::applyJacobianInverse(
					  Teuchos::ParameterList& params, 
					  const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  // Convert input, result to multivectors
  Teuchos::RCP<NOX::Abstract::MultiVector> mv_input = 
    input.createMultiVector(1, NOX::DeepCopy);
  Teuchos::RCP<NOX::Abstract::MultiVector> mv_result = 
    result.createMultiVector(1, NOX::DeepCopy);

  // Call multivector version of applyJacobianInverse
  NOX::Abstract::Group::ReturnType status = 
    applyJacobianInverseMultiVector(params, *mv_input, *mv_result);

  // Copy result
  result = (*mv_result)[0];

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::applyJacobianMultiVector(
				     const NOX::Abstract::MultiVector& input,
				     NOX::Abstract::MultiVector& result) const 
{
  std::string callingFunction = 
    "LOCA::Hopf::MooreSpence::ExtendedGroup::applyJacobianMultiVector()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  if (!isJacobian()) {
    globalData->locaErrorCheck->throwError(callingFunction,
					   "Called with invalid Jacobian!");
  }

  // Cast vectors to Hopf vectors
  const LOCA::Hopf::MooreSpence::ExtendedMultiVector& hopf_input = 
    dynamic_cast<const LOCA::Hopf::MooreSpence::ExtendedMultiVector&>(input);
  LOCA::Hopf::MooreSpence::ExtendedMultiVector& hopf_result = 
    dynamic_cast<LOCA::Hopf::MooreSpence::ExtendedMultiVector&>(result);

  // Get constant references to input vector components
  Teuchos::RCP<const NOX::Abstract::MultiVector> input_x = 
    hopf_input.getXMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector> input_y = 
    hopf_input.getRealEigenMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector> input_z = 
    hopf_input.getImagEigenMultiVec();
  Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> input_w =
    hopf_input.getFrequencies();
  Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> input_p =
    hopf_input.getBifParams();

  // Get non-constant references to result vector components
  Teuchos::RCP<NOX::Abstract::MultiVector> result_x = 
    hopf_result.getXMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector> result_y = 
    hopf_result.getRealEigenMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector> result_z = 
    hopf_result.getImagEigenMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> result_w = 
    hopf_result.getFrequencies();
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> result_p = 
    hopf_result.getBifParams();

  // Temporary vectors
  Teuchos::RCP<NOX::Abstract::MultiVector> tmp_real = 
    input_y->clone(NOX::ShapeCopy);
  Teuchos::RCP<NOX::Abstract::MultiVector> tmp_imag = 
    input_z->clone(NOX::ShapeCopy);

  // verify underlying Jacobian is valid
  if (!grpPtr->isJacobian()) {
    status = grpPtr->computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  // compute J*X
  status = grpPtr->applyJacobianMultiVector(*input_x, *result_x);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							   finalStatus,
							   callingFunction);

  // compute J*X + P*dR/dp
  result_x->update(Teuchos::NO_TRANS, 1.0, *(dfdpMultiVec->getXMultiVec()), 
		   *input_p);

  // compute (J+iwB)*(Y+iZ)
  status = grpPtr->applyComplexMultiVector(*input_y, *input_z, *result_y,
					   *result_z);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  // compute (J+iwB)*(Y+iZ) + P*d(J+iwB)*(y+iz)/dp
  result_y->update(Teuchos::NO_TRANS, 1.0, 
		   *(dfdpMultiVec->getRealEigenMultiVec()), *input_p);
  result_z->update(Teuchos::NO_TRANS, 1.0, 
		   *(dfdpMultiVec->getImagEigenMultiVec()), *input_p);

  // compute (J+iwB)*(Y+iZ) + P*d(J+iwB)*(y+iz)/dp + iW*B*(y+iz)
  result_y->update(Teuchos::NO_TRANS, 1.0, *minusMassTimesZ, *input_w);
  result_z->update(Teuchos::NO_TRANS, 1.0, *massTimesY, *input_w);

  // compute (d(J+iwB)(y+iz)/dx)*X
  status = grpPtr->computeDCeDxa(*(xVec->getRealEigenVec()), 
				 *(xVec->getImagEigenVec()), 
				 xVec->getFrequency(),
				 *input_x,
				 *(fVec->getRealEigenVec()),
				 *(fVec->getImagEigenVec()),
				 *tmp_real, *tmp_imag);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  
  // (d(J+iwB)*(y+iz)/dx)*X+(J+iwB)*(Y+iZ) P*d(J+iwB)*(y+iz)/dp+iW*B*(y+iz)
  result_y->update(1.0, *tmp_real, 1.0);
  result_z->update(1.0, *tmp_imag, 1.0);

  // compute phi^T*Y, phi^T*Z
  lTransNorm(*input_y, *result_w);
  lTransNorm(*input_z, *result_p);

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::applyJacobianTransposeMultiVector(
				     const NOX::Abstract::MultiVector& input,
				     NOX::Abstract::MultiVector& result) const 
{
  globalData->locaErrorCheck->throwError(
		  "LOCA::Hopf::MooreSpence::ExtendedGroup::applyJacobianTransposeMultiVector()",
		  "Method not implemented!");

  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::applyJacobianInverseMultiVector( 
				     Teuchos::ParameterList& params, 	
				     const NOX::Abstract::MultiVector& input,
				     NOX::Abstract::MultiVector& result) const 
{
  // Cast vectors to Hopf vectors
  const LOCA::Hopf::MooreSpence::ExtendedMultiVector& hopf_input = 
    dynamic_cast<const LOCA::Hopf::MooreSpence::ExtendedMultiVector&>(input);
  LOCA::Hopf::MooreSpence::ExtendedMultiVector& hopf_result = 
    dynamic_cast<LOCA::Hopf::MooreSpence::ExtendedMultiVector&>(result);
  
  return solverStrategy->solve(params, hopf_input, hopf_result);
}

bool
LOCA::Hopf::MooreSpence::ExtendedGroup::isF() const 
{
  return isValidF;
}

bool
LOCA::Hopf::MooreSpence::ExtendedGroup::isJacobian() const 
{
  return isValidJacobian;
}

bool
LOCA::Hopf::MooreSpence::ExtendedGroup::isGradient() const 
{
  return false;
}

bool
LOCA::Hopf::MooreSpence::ExtendedGroup::isNewton() const 
{
  return isValidNewton;
}
  
const NOX::Abstract::Vector&
LOCA::Hopf::MooreSpence::ExtendedGroup::getX() const 
{
  return *xVec;
}

const NOX::Abstract::Vector&
LOCA::Hopf::MooreSpence::ExtendedGroup::getF() const 
{
  return *fVec;
}

double
LOCA::Hopf::MooreSpence::ExtendedGroup::getNormF() const 
{
  return fVec->norm();
}

const NOX::Abstract::Vector&
LOCA::Hopf::MooreSpence::ExtendedGroup::getGradient() const 
{
  globalData->locaErrorCheck->throwError(
	      "LOCA::Hopf::MooreSpence::ExtendedGroup::getGradient()",
	      " - not implemented");
  return getNewton();
}

const NOX::Abstract::Vector&
LOCA::Hopf::MooreSpence::ExtendedGroup::getNewton() const 
{
  return *newtonVec;
}
  
Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Hopf::MooreSpence::ExtendedGroup::getXPtr() const 
{
  return xVec;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Hopf::MooreSpence::ExtendedGroup::getFPtr() const 
{
  return fVec;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Hopf::MooreSpence::ExtendedGroup::getGradientPtr() const 
{
  globalData->locaErrorCheck->throwError(
	      "LOCA::Hopf::MooreSpence::ExtendedGroup::getGradientPtr()",
	      " - not implemented");
  return getNewtonPtr();
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Hopf::MooreSpence::ExtendedGroup::getNewtonPtr() const 
{
  return newtonVec;
}

double
LOCA::Hopf::MooreSpence::ExtendedGroup::getNormNewtonSolveResidual() const 
{
  std::string callingFunction = 
    "LOCA::Hopf::MooreSpence::ExtendedGroup::getNormNewtonSolveResidual()";
  NOX::Abstract::Group::ReturnType finalStatus;
  LOCA::Hopf::MooreSpence::ExtendedVector residual = *fVec;
  
  finalStatus = applyJacobian(*newtonVec, residual);
  globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

  residual.update(1.0, *fVec, 1.0);
  return residual.norm();
}

Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup>
LOCA::Hopf::MooreSpence::ExtendedGroup::getUnderlyingGroup() const
{
  return grpPtr;
}

Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
LOCA::Hopf::MooreSpence::ExtendedGroup::getUnderlyingGroup()
{
  return grpPtr;
}

void
LOCA::Hopf::MooreSpence::ExtendedGroup::copy(
					    const NOX::Abstract::Group& src)
{
  const LOCA::Hopf::MooreSpence::ExtendedGroup& source = 
    dynamic_cast<const LOCA::Hopf::MooreSpence::ExtendedGroup&>(src);

  // Protect against A = A
  if (this != &source) {
    
    // Copy values
    globalData = source.globalData;
    parsedParams = source.parsedParams;
    hopfParams = source.hopfParams;
    grpPtr->copy(*(source.grpPtr));
    xMultiVec = source.xMultiVec;
    fMultiVec = source.fMultiVec;
    newtonMultiVec = source.newtonMultiVec;
    *lengthMultiVec = *source.lengthMultiVec;
    *massTimesY = *source.massTimesY;
    *minusMassTimesZ = *source.minusMassTimesZ;
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
      globalData->locaFactory->createMooreSpenceHopfSolverStrategy(
								 parsedParams,
								 hopfParams);
  }
}

void
LOCA::Hopf::MooreSpence::ExtendedGroup::setParamsMulti(
			  const std::vector<int>& paramIDs, 
			  const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  grpPtr->setParamsMulti(paramIDs, vals);
  for (unsigned int i=0; i<paramIDs.size(); i++) {
    if (paramIDs[i] == bifParamID[0])
      setBifParam(vals(i,0));
  }
}

void
LOCA::Hopf::MooreSpence::ExtendedGroup::setParams(
					      const LOCA::ParameterVector& p) 
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  grpPtr->setParams(p);
  setBifParam(p[bifParamID[0]]);
}

void
LOCA::Hopf::MooreSpence::ExtendedGroup::setParam(int paramID, double val)
{
  if (paramID == bifParamID[0])
    setBifParam(val);
  else
    grpPtr->setParam(paramID, val);
}

void
LOCA::Hopf::MooreSpence::ExtendedGroup::setParam(std::string paramID, double val)
{
  const LOCA::ParameterVector& pVec = grpPtr->getParams();
  if (pVec.getIndex(paramID) == bifParamID[0])
    setBifParam(val);
  else
    grpPtr->setParam(paramID, val);
}

const LOCA::ParameterVector&
LOCA::Hopf::MooreSpence::ExtendedGroup::getParams() const 
{
  return grpPtr->getParams();
}

double
LOCA::Hopf::MooreSpence::ExtendedGroup::getParam(int paramID) const
{
  return grpPtr->getParam(paramID);
}

double
LOCA::Hopf::MooreSpence::ExtendedGroup::getParam(std::string paramID) const
{
  return grpPtr->getParam(paramID);
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::computeDfDpMulti(
					    const std::vector<int>& paramIDs, 
					    NOX::Abstract::MultiVector& dfdp, 
					    bool isValid_F)
{
   std::string callingFunction = 
    "LOCA::Hopf::MooreSpence::ExtendedGroup::computeDfDpMulti()";
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Cast dfdp to Hopf vector
  LOCA::Hopf::MooreSpence::ExtendedMultiVector& hopf_dfdp = 
    dynamic_cast<LOCA::Hopf::MooreSpence::ExtendedMultiVector&>(dfdp);

  // Compute df/dp
  status = grpPtr->computeDfDpMulti(paramIDs, *hopf_dfdp.getXMultiVec(),
				    isValid_F);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  // Compute d(J+iwB)*(y+iz)/dp
  status = grpPtr->computeDCeDp(paramIDs, 
				*(xVec->getRealEigenVec()),
				*(xVec->getImagEigenVec()),
				xVec->getFrequency(),
				*(hopf_dfdp.getRealEigenMultiVec()), 
				*(hopf_dfdp.getImagEigenMultiVec()),
				isValid_F);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  // Set frequency & parameter components
  if (!isValid_F) {
    hopf_dfdp.getScalar(0,0) = lTransNorm(*(xVec->getRealEigenVec()));
    hopf_dfdp.getScalar(1,0) = lTransNorm(*(xVec->getImagEigenVec()));
  }
  for (int i=0; i<dfdp.numVectors()-1; i++) {
    hopf_dfdp.getScalar(0,i+1) = 0.0;
    hopf_dfdp.getScalar(1,i+1) = 0.0;
  }

  return finalStatus;
}

void
LOCA::Hopf::MooreSpence::ExtendedGroup::preProcessContinuationStep(
			 LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  grpPtr->preProcessContinuationStep(stepStatus);
}

void
LOCA::Hopf::MooreSpence::ExtendedGroup::postProcessContinuationStep(
			 LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  grpPtr->postProcessContinuationStep(stepStatus);
}

void
LOCA::Hopf::MooreSpence::ExtendedGroup::projectToDraw(
					       const NOX::Abstract::Vector& x,
					       double *px) const
{
  const LOCA::Hopf::MooreSpence::ExtendedVector& mx = 
    dynamic_cast<const LOCA::Hopf::MooreSpence::ExtendedVector&>(x);

  grpPtr->projectToDraw(*(mx.getXVec()), px);
  px[grpPtr->projectToDrawDimension()] = mx.getFrequency();
  px[grpPtr->projectToDrawDimension()+1] = mx.getBifParam();
}

int
LOCA::Hopf::MooreSpence::ExtendedGroup::projectToDrawDimension() const
{
  return grpPtr->projectToDrawDimension() + 2;
}

void
LOCA::Hopf::MooreSpence::ExtendedGroup::printSolution(
						  const double conParam) const 
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() << 
      "LOCA::Hopf::MooreSpence::ExtendedGroup::printSolution\n";

    globalData->locaUtils->out() << "Hopf Point located at: " << 
      globalData->locaUtils->sciformat(conParam) << "   " << 
      globalData->locaUtils->sciformat(xVec->getBifParam()) << "   " << 
      globalData->locaUtils->sciformat(xVec->getFrequency()) << std::endl;

    globalData->locaUtils->out() << 
      "\tPrinting Solution Vector for conParam = " << 
      globalData->locaUtils->sciformat(conParam) << std::endl;
  }
  grpPtr->printSolution(conParam);
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() << 
      "\tPrinting Real Component of Eigenvector for bif param = " << 
      globalData->locaUtils->sciformat(xVec->getBifParam()) << std::endl;
  }
  grpPtr->printSolution(*(xVec->getRealEigenVec()), xVec->getBifParam());
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() << 
      "\tPrinting Imaginary Component of Eigenvector for frequency = " << 
      globalData->locaUtils->sciformat(xVec->getFrequency()) << std::endl;
  }
  grpPtr->printSolution(*(xVec->getImagEigenVec()), xVec->getFrequency());
}

void
LOCA::Hopf::MooreSpence::ExtendedGroup::printSolution(
					     const NOX::Abstract::Vector& x_,
					     const double conParam) const 
{
  const LOCA::Hopf::MooreSpence::ExtendedVector& hopf_x = 
    dynamic_cast<const LOCA::Hopf::MooreSpence::ExtendedVector&>(x_);

  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() << 
      "LOCA::Hopf::MooreSpence::ExtendedGroup::printSolution\n";

    globalData->locaUtils->out() << "Hopf Point located at: " << 
      globalData->locaUtils->sciformat(conParam) << "   " << 
      globalData->locaUtils->sciformat(hopf_x.getBifParam()) << "   " << 
      globalData->locaUtils->sciformat(hopf_x.getFrequency())<< std::endl;

     globalData->locaUtils->out() << 
       "\tPrinting Solution Vector for conParam = " << 
       globalData->locaUtils->sciformat(conParam) << std::endl;
  }
  grpPtr->printSolution(*hopf_x.getXVec(), conParam);
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() << 
      "\tPrinting Real Component of Eigenvector for bif param = " << 
      globalData->locaUtils->sciformat(hopf_x.getBifParam()) << std::endl;
  }
  grpPtr->printSolution(*hopf_x.getRealEigenVec(), hopf_x.getBifParam());
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() << 
      "\tPrinting Imaginary Component of Eigenvector for frequency = " << 
      globalData->locaUtils->sciformat(hopf_x.getFrequency()) << std::endl;
  }
  grpPtr->printSolution(*hopf_x.getImagEigenVec(), hopf_x.getFrequency());
}

double
LOCA::Hopf::MooreSpence::ExtendedGroup::getBifParam() const 
{
  return grpPtr->getParam(bifParamID[0]);
}

double
LOCA::Hopf::MooreSpence::ExtendedGroup::getFrequency() const 
{
  return xVec->getFrequency();
}

double
LOCA::Hopf::MooreSpence::ExtendedGroup::lTransNorm(
					const NOX::Abstract::Vector& n) const
{
  return lengthVec->innerProduct(n) / lengthVec->length();
}

void
LOCA::Hopf::MooreSpence::ExtendedGroup::lTransNorm(
			const NOX::Abstract::MultiVector& n,
			NOX::Abstract::MultiVector::DenseMatrix& result) const
{
  n.multiply(1.0 / lengthVec->length(), *lengthMultiVec, result);
}

void
LOCA::Hopf::MooreSpence::ExtendedGroup::setBifParam(double param) 
{
  grpPtr->setParam(bifParamID[0], param);
  xVec->getBifParam() = param;

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

void
LOCA::Hopf::MooreSpence::ExtendedGroup::setupViews()
{
  index_f[0] = 0;
  index_dfdp[0] = 1;
  
  xVec = xMultiVec.getColumn(0);
  fVec = fMultiVec.getColumn(0);
  newtonVec = newtonMultiVec.getColumn(0);
  lengthVec = Teuchos::rcp(&(*lengthMultiVec)[0],false);

  ffMultiVec = Teuchos::rcp_dynamic_cast<LOCA::Hopf::MooreSpence::ExtendedMultiVector>(fMultiVec.subView(index_f),true);

  dfdpMultiVec = Teuchos::rcp_dynamic_cast<LOCA::Hopf::MooreSpence::ExtendedMultiVector>(fMultiVec.subView(index_dfdp),true);

}

void
LOCA::Hopf::MooreSpence::ExtendedGroup::init(bool perturbSoln, 
					     double perturbSize)
{
  xVec->getBifParam() = getBifParam();

   // Rescale and rotate complex eigenvector so normalization condition is met
  double ldy = lTransNorm(*(xVec->getRealEigenVec()));
  double ldz = lTransNorm(*(xVec->getImagEigenVec()));

  if (fabs(ldy) < 1.0e-8) {
    globalData->locaErrorCheck->throwError(
		   "LOCA::Hopf::MooreSpence::ExtendedGroup::init()",
		   "Real component of eigenvector cannot be orthogonal to length-scaling vector ");
  }

  double denom = ldy*ldy + ldz*ldz;
  double a =  ldy/denom;
  double b = -ldz/denom;
  Teuchos::RCP<NOX::Abstract::Vector> y_tmp = 
    xVec->getRealEigenVec()->clone();

  // y <- a*y - b*z
  xVec->getRealEigenVec()->update(-b, *(xVec->getImagEigenVec()), a);

  // z <- a*z + b*y
  xVec->getImagEigenVec()->update(b, *y_tmp, a);

  if (perturbSoln) {
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
     globalData->locaUtils->out() << 
       "\tIn LOCA::Hopf::MooreSpence::ExtendedGroup::init(), " << 
       "applying random perturbation to initial solution of size: " << 
       globalData->locaUtils->sciformat(perturbSize) << std::endl;
    }
    Teuchos::RCP<NOX::Abstract::Vector> perturb = 
      xVec->getXVec()->clone(NOX::ShapeCopy);
    perturb->random();
    perturb->scale(*(xVec->getXVec()));
    xVec->getXVec()->update(perturbSize, *perturb, 1.0);
    grpPtr->setX(*(xVec->getXVec()));
  }
}

