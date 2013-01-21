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

#include "LOCA_TurningPoint_MinimallyAugmented_ModifiedConstraint.H"
#include "LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H"
#include "LOCA_BorderedSolver_AbstractStrategy.H"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Factory.H"
#include "NOX_Utils.H"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_BorderedSolver_JacobianOperator.H"

LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::
ModifiedConstraint(
    const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
    const Teuchos::RCP<Teuchos::ParameterList>& tpParams,
    const Teuchos::RCP<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup>& g,
    int bif_param) :
  Constraint(global_data, topParams, tpParams, g, bif_param),
  w_vector_update(a_vector->clone(NOX::ShapeCopy)),
  v_vector_update(a_vector->clone(NOX::ShapeCopy)),
  w_residual(a_vector->clone(NOX::ShapeCopy)),
  v_residual(a_vector->clone(NOX::ShapeCopy)),
  deltaX(a_vector->clone(NOX::ShapeCopy)),
  sigma1(1, 1),
  sigma2(1, 1),
  deltaP(0),
  isFirstSolve(true),
  includeNewtonTerms(false)
{
  // zero out null vector updates
  w_vector_update->init(0.0);
  v_vector_update->init(0.0);

  includeNewtonTerms = tpParams->get("Include Newton Terms", false);
}

LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::
ModifiedConstraint(
     const LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint& source, 
     NOX::CopyType type) : 
  Constraint(source, type),
  w_vector_update(source.w_vector_update->clone(type)),
  v_vector_update(source.v_vector_update->clone(type)),
  w_residual(source.w_residual->clone(type)),
  v_residual(source.v_residual->clone(type)),
  deltaX(source.deltaX->clone(type)),
  sigma1(source.sigma1),
  sigma2(source.sigma2),
  deltaP(source.deltaP),
  isFirstSolve(source.isFirstSolve),
  includeNewtonTerms(source.includeNewtonTerms)
{
}

LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::
~ModifiedConstraint()
{
}

void
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::
copy(const LOCA::MultiContinuation::ConstraintInterface& src)
{
  const LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint& source = 
  dynamic_cast<const LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint&>(src);

  if (this != &source) {
    Constraint::copy(source);

    *w_vector_update = *(source.w_vector_update);
    *v_vector_update = *(source.v_vector_update);
    *w_residual = *(source.w_residual);
    *v_residual = *(source.v_residual);
    *deltaX = *(source.deltaX);
    sigma1.assign(source.sigma1);
    sigma2.assign(source.sigma2);
    deltaP = source.deltaP;
    isFirstSolve = source.isFirstSolve;
    includeNewtonTerms = source.includeNewtonTerms;
  }
}

Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::
clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new ModifiedConstraint(*this, type));
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::
computeConstraints()
{
  if (isValidConstraints)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction = 
    "LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::computeConstraints()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // Compute J
  status = grpPtr->computeJacobian();
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							   finalStatus,
							   callingFunction);

  // Set up bordered systems
  Teuchos::RCP<const LOCA::BorderedSolver::JacobianOperator> op =
    Teuchos::rcp(new  LOCA::BorderedSolver::JacobianOperator(grpPtr));
  borderedSolver->setMatrixBlocksMultiVecConstraint(op, 
						    a_vector, 
						    b_vector, 
						    Teuchos::null);

  // Get linear solver parameters
  Teuchos::RCP<Teuchos::ParameterList> linear_solver_params =
    parsedParams->getSublist("Linear Solver");

  // Solve for w and v
  if (isFirstSolve) {

    std::cout << "solving for base w,v..." << std::endl;

    // Create RHS
    NOX::Abstract::MultiVector::DenseMatrix one(1,1);
    if (nullVecScaling == NVS_OrderN)
      one(0,0) = dn;
    else
      one(0,0) = 1.0;

    // Compute sigma_1 and right null vector v
    status = borderedSolver->initForSolve();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
    status = borderedSolver->applyInverse(*linear_solver_params, 
					  NULL, 
					  &one, 
					  *v_vector, 
					  sigma1);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);

    // Compute sigma_2 and left null vector w
    if (!isSymmetric) {
      status = borderedSolver->initForTransposeSolve();
      finalStatus = 
	globalData->locaErrorCheck->combineAndCheckReturnTypes(
							    status, 
							    finalStatus,
							    callingFunction);
      status = borderedSolver->applyInverseTranspose(*linear_solver_params, 
						     NULL, 
						     &one, 
						     *w_vector, 
						     sigma2);
      finalStatus = 
	globalData->locaErrorCheck->combineAndCheckReturnTypes(
							    status, 
							    finalStatus,
							    callingFunction);

    }
    else {
      *w_vector = *v_vector;
      sigma2.assign(sigma1);
    }

    isFirstSolve = false;
  }

  // solve for updates to w and v
  else {

    std::cout << "solving for updates..." << std::endl;

    // Compute J*v + a*sigma_1
    status = grpPtr->applyJacobianMultiVector(*v_vector, *v_residual);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
    v_residual->update(Teuchos::NO_TRANS, 1.0, *a_vector, sigma1, 0.0);
    
    // Compute b^T*v - n
    NOX::Abstract::MultiVector::DenseMatrix sigma1_residual(1,1);
    v_vector->multiply(1.0, *b_vector, sigma1_residual);
    if (nullVecScaling == NVS_OrderN)
      sigma1_residual(0,0) -= dn;
    else
      sigma1_residual(0,0) -= 1.0;

    if (includeNewtonTerms) {

      // Compute (Jv)_x*dx
      Teuchos::RCP<NOX::Abstract::MultiVector> Jv_x_dx = 
	deltaX->clone(NOX::ShapeCopy);
      status = grpPtr->computeDJnDxaMulti((*v_vector)[0], *deltaX, *Jv_x_dx);
      finalStatus = 
	globalData->locaErrorCheck->combineAndCheckReturnTypes(
							     status, 
							     finalStatus,
							     callingFunction);

      // Compute (Jv)_p
      Teuchos::RCP<NOX::Abstract::MultiVector> Jv_p1 = 
	deltaX->clone(2);
      std::vector<int> idx(1); idx[0] = 0;
      Teuchos::RCP<NOX::Abstract::MultiVector> Jv_p = 
	Jv_p1->subView(idx);
      status = grpPtr->computeDJnDpMulti(bifParamID, (*v_vector)[0], *Jv_p1, 
					 false);
      finalStatus = 
	globalData->locaErrorCheck->combineAndCheckReturnTypes(
							     status, 
							     finalStatus,
							     callingFunction);
      // compute v_residual += (Jv)_x*dx + (Jv)_p*dp
      v_residual->update(1.0, *Jv_x_dx, deltaP, *Jv_p, 1.0);
      
      // Compute J
      status = grpPtr->computeJacobian();
      finalStatus = 
	globalData->locaErrorCheck->combineAndCheckReturnTypes(
							     status, 
							     finalStatus,
							     callingFunction);
    }

    // Compute update to sigma_1 and right null vector v
    NOX::Abstract::MultiVector::DenseMatrix sigma1_update(1,1);
    status = borderedSolver->initForSolve();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(
							   status, 
							   finalStatus,
							   callingFunction);
    status = borderedSolver->applyInverse(*linear_solver_params, 
					  v_residual.get(), 
					  &sigma1_residual, 
					  *v_vector_update, 
					  sigma1_update);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(
							   status, 
							   finalStatus,
							   callingFunction);

    // apply updates
    v_vector->update(-1.0, *v_vector_update, 1.0);
    sigma1(0,0) -= sigma1_update(0,0);
    
    if (!isSymmetric) {
      // Compute J^T*w + b*sigma_w
      status = grpPtr->applyJacobianTransposeMultiVector(*w_vector, 
							 *w_residual);
      finalStatus = 
	globalData->locaErrorCheck->combineAndCheckReturnTypes(
							     status, 
							     finalStatus,
							     callingFunction);
      w_residual->update(Teuchos::NO_TRANS, 1.0, *b_vector, sigma2, 0.0);

      // Compute a^T*w - n
      NOX::Abstract::MultiVector::DenseMatrix sigma2_residual(1,1);
      w_vector->multiply(1.0, *a_vector, sigma2_residual);
      if (nullVecScaling == NVS_OrderN)
	sigma2_residual(0,0) -= dn;
      else
	sigma2_residual(0,0) -= 1.0;

      if (includeNewtonTerms) {

	// Compute (J^T*w)_x*dx
	Teuchos::RCP<NOX::Abstract::MultiVector> Jtw_x_dx = 
	  deltaX->clone(NOX::ShapeCopy);
	status = grpPtr->computeDwtJnDx((*w_vector)[0], (*deltaX)[0], 
					(*Jtw_x_dx)[0]);
	finalStatus = 
	  globalData->locaErrorCheck->combineAndCheckReturnTypes(
							     status, 
							     finalStatus,
							     callingFunction);

	// Compute (J^T*w)_p
	Teuchos::RCP<NOX::Abstract::MultiVector> Jtw_p1 = 
	  deltaX->clone(2);
	std::vector<int> idx(1); idx[0] = 0;
	Teuchos::RCP<NOX::Abstract::MultiVector> Jtw_p = 
	  Jtw_p1->subView(idx);
	status = grpPtr->computeDwtJDp(bifParamID, (*w_vector)[0], *Jtw_p1, 
				       false);
	finalStatus = 
	  globalData->locaErrorCheck->combineAndCheckReturnTypes(
							     status, 
							     finalStatus,
							     callingFunction);
	// compute w_residual += (J^T*w)_x*dx + (J^T*w)_p*dp
	w_residual->update(1.0, *Jtw_x_dx, deltaP, *Jtw_p, 1.0);
	
	// Compute J
	status = grpPtr->computeJacobian();
	finalStatus = 
	  globalData->locaErrorCheck->combineAndCheckReturnTypes(
							    status, 
							    finalStatus,
							    callingFunction);
      }

      // Compute update to sigma_2 and left null vector w
      NOX::Abstract::MultiVector::DenseMatrix sigma2_update(1,1);
      status = borderedSolver->initForTransposeSolve();
      finalStatus = 
	globalData->locaErrorCheck->combineAndCheckReturnTypes(
							     status, 
							     finalStatus,
							     callingFunction);
      status = borderedSolver->applyInverseTranspose(*linear_solver_params, 
						     w_residual.get(), 
						     &sigma2_residual, 
						     *w_vector_update, 
						     sigma2_update);
      finalStatus = 
	globalData->locaErrorCheck->combineAndCheckReturnTypes(
							     status, 
							     finalStatus,
							     callingFunction);

      // apply updates
      w_vector->update(-1.0, *w_vector_update, 1.0);
      sigma2(0,0) -= sigma2_update(0,0);

    }
    else {
      *w_vector = *v_vector;
      sigma2.assign(sigma1);
    }
    
  }
  
  // Compute sigma = -w^T*J*v
  status = grpPtr->applyJacobianMultiVector(*v_vector, *Jv_vector);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							   finalStatus,
							   callingFunction);
  if (!isSymmetric) {
    status = grpPtr->applyJacobianTransposeMultiVector(*w_vector, *Jtw_vector);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }
  else
    *Jtw_vector = *Jv_vector;
  Jv_vector->multiply(-1.0, *w_vector, constraints);

  // Scale sigma
  double w_norm = (*w_vector)[0].norm();
  double v_norm = (*v_vector)[0].norm();
  double Jv_norm = (*Jv_vector)[0].norm();
  double Jtw_norm = (*Jtw_vector)[0].norm();
  if (nullVecScaling == NVS_OrderN)
    sigma_scale = dn;
  else
    sigma_scale = 1.0;
  constraints.scale(1.0/sigma_scale);

  if (globalData->locaUtils->isPrintType(NOX::Utils::OuterIteration)) {
    globalData->locaUtils->out() << "\n\t||Right null vector v|| = " 
				 << globalData->locaUtils->sciformat(v_norm);
    globalData->locaUtils->out() << "\n\t||Left null vector w|| = " 
				 << globalData->locaUtils->sciformat(w_norm);
    globalData->locaUtils->out() << "\n\t||Jv|| = " 
				 << globalData->locaUtils->sciformat(Jv_norm);
    globalData->locaUtils->out() << "\n\t||J^T*w|| = " 
				 << globalData->locaUtils->sciformat(Jtw_norm);
    globalData->locaUtils->out() << 
      "\n\tRight estimate for singularity of Jacobian (sigma1) = " << 
      globalData->locaUtils->sciformat(sigma1(0,0));
    globalData->locaUtils->out() << 
      "\n\tLeft estimate for singularity of Jacobian (sigma2) = " << 
      globalData->locaUtils->sciformat(sigma2(0,0));
    globalData->locaUtils->out() << 
      "\n\tFinal Estimate for singularity of Jacobian (sigma) = " << 
      globalData->locaUtils->sciformat(constraints(0,0)) << std::endl;
  }

  isValidConstraints = true;

  // Update a and b if requested
  if (updateVectorsEveryIteration) {
    if (globalData->locaUtils->isPrintType(NOX::Utils::OuterIteration)) {
      globalData->locaUtils->out() << 
	"\n\tUpdating null vectors for the next nonlinear iteration" << 
	std::endl;
    }
    *a_vector = *w_vector;
    *b_vector = *v_vector;

    scaleNullVectors((*a_vector)[0],(*b_vector)[0]);
  }

  return finalStatus;
}

void
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::
preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  Constraint::preProcessContinuationStep(stepStatus);

    // zero out null vector updates
    w_vector_update->init(0.0);
    v_vector_update->init(0.0);

    isFirstSolve = true;
}

void
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::
postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  Constraint::postProcessContinuationStep(stepStatus);

  if (stepStatus == LOCA::Abstract::Iterator::Successful) {

    // zero out null vector updates
    w_vector_update->init(0.0);
    v_vector_update->init(0.0);

    isFirstSolve = true;
  }
}

void
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::
setNewtonUpdates(const NOX::Abstract::Vector& dx, double dp, double step)
{
  (*deltaX)[0].update(step, dx, 0.0);
  deltaP = step*dp;
}


