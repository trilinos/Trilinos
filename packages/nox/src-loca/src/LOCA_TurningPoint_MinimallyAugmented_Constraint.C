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

#include "LOCA_TurningPoint_MinimallyAugmented_Constraint.H"
#include "LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H"
#include "LOCA_BorderedSolver_AbstractStrategy.H"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Factory.H"
#include "NOX_Utils.H"
#include "Teuchos_ParameterList.hpp"

LOCA::TurningPoint::MinimallyAugmented::Constraint::
Constraint(
    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
    const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
    const Teuchos::RefCountPtr<Teuchos::ParameterList>& tpParams,
    const Teuchos::RefCountPtr<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup>& g,
    bool is_symmetric,
    const NOX::Abstract::Vector& a,
    const NOX::Abstract::Vector* b,
    int bif_param) :
  globalData(global_data),
  parsedParams(topParams),
  turningPointParams(tpParams),
  grpPtr(g),
  a_vector(a.createMultiVector(1, NOX::DeepCopy)),
  b_vector(),
  w_vector(a.createMultiVector(1, NOX::ShapeCopy)),
  v_vector(a.createMultiVector(1, NOX::ShapeCopy)),
  Jv_vector(a.createMultiVector(1, NOX::ShapeCopy)),
  sigma_x(a.createMultiVector(1, NOX::ShapeCopy)),
  constraints(1, 1),
  borderedSolver(),
  dn(static_cast<double>(a_vector->length())),
  sigma_scale(1.0),
  isSymmetric(is_symmetric),
  isValidConstraints(false),
  isValidDX(false),
  bifParamID(1),
  updateVectorsEveryContinuationStep(true),
  updateVectorsEveryIteration(false)
{
  // Instantiate bordered solvers
  borderedSolver = 
    globalData->locaFactory->createBorderedSolverStrategy(parsedParams,
							  turningPointParams);
  if (!isSymmetric) {
    b_vector = b->createMultiVector(1, NOX::DeepCopy);
  }
  else {
    b_vector = a_vector->clone(NOX::DeepCopy);
  }

  // Options
  updateVectorsEveryContinuationStep = 
    turningPointParams->get(
			       "Update Null Vectors Every Continuation Step", 
			       true);
  updateVectorsEveryIteration = 
    turningPointParams->get(
			      "Update Null Vectors Every Nonlinear Iteration", 
			      false);
}

LOCA::TurningPoint::MinimallyAugmented::Constraint::
Constraint(const LOCA::TurningPoint::MinimallyAugmented::Constraint& source, 
	   NOX::CopyType type) : 
  globalData(source.globalData),
  parsedParams(source.parsedParams),
  turningPointParams(source.turningPointParams),
  grpPtr(Teuchos::null),
  a_vector(source.a_vector->clone(type)),
  b_vector(source.b_vector->clone(type)),
  w_vector(source.w_vector->clone(type)),
  v_vector(source.v_vector->clone(type)),
  Jv_vector(source.Jv_vector->clone(type)),
  sigma_x(source.sigma_x->clone(type)),
  constraints(source.constraints),
  borderedSolver(),
  dn(source.dn),
  sigma_scale(source.sigma_scale),
  isSymmetric(source.isSymmetric),
  isValidConstraints(false),
  isValidDX(false),
  bifParamID(source.bifParamID),
  updateVectorsEveryContinuationStep(source.updateVectorsEveryContinuationStep),
  updateVectorsEveryIteration(source.updateVectorsEveryIteration)
{
  if (source.isValidConstraints && type == NOX::DeepCopy)
    isValidConstraints = true;
  if (source.isValidDX && type == NOX::DeepCopy)
    isValidDX = true;

  // Instantiate bordered solvers
  borderedSolver = 
    globalData->locaFactory->createBorderedSolverStrategy(parsedParams,
							  turningPointParams);

  // We don't explicitly copy the group because the constrained group
  // will do that
}

LOCA::TurningPoint::MinimallyAugmented::Constraint::
~Constraint()
{
}

void
LOCA::TurningPoint::MinimallyAugmented::Constraint::
setGroup(const Teuchos::RefCountPtr<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup>& g)
{
  grpPtr = g;
}

Teuchos::RefCountPtr<const NOX::Abstract::Vector>
LOCA::TurningPoint::MinimallyAugmented::Constraint::
getLeftNullVec() const
{
  return Teuchos::rcp(&(*w_vector)[0], false);
}

Teuchos::RefCountPtr<const NOX::Abstract::Vector>
LOCA::TurningPoint::MinimallyAugmented::Constraint::
getRightNullVec() const
{
  return Teuchos::rcp(&(*v_vector)[0], false);
}

double
LOCA::TurningPoint::MinimallyAugmented::Constraint::
getSigma() const
{
  return constraints(0,0);
}

void
LOCA::TurningPoint::MinimallyAugmented::Constraint::
copy(const LOCA::MultiContinuation::ConstraintInterface& src)
{
  const LOCA::TurningPoint::MinimallyAugmented::Constraint& source = 
  dynamic_cast<const LOCA::TurningPoint::MinimallyAugmented::Constraint&>(src);

  if (this != &source) {
    globalData = source.globalData;
    parsedParams = source.parsedParams;
    turningPointParams = source.turningPointParams;
    *a_vector = *(source.a_vector);
    *b_vector = *(source.b_vector);
    *w_vector = *(source.w_vector);
    *v_vector = *(source.v_vector);
    *Jv_vector = *(source.Jv_vector);
    *sigma_x = *(source.sigma_x);
    constraints.assign(source.constraints);
    dn = source.dn;
    sigma_scale = source.sigma_scale;
    isSymmetric = source.isSymmetric;
    isValidConstraints = source.isValidConstraints;
    isValidDX = source.isValidDX;
    bifParamID = source.bifParamID;
    updateVectorsEveryContinuationStep = 
      source.updateVectorsEveryContinuationStep;
    updateVectorsEveryIteration = 
      source.updateVectorsEveryIteration;

    // Instantiate bordered solvers
    borderedSolver = 
      globalData->locaFactory->createBorderedSolverStrategy(
							 parsedParams,
							 turningPointParams);

    // We don't explicitly copy the group because the constrained group
    // will do that
  }
}

Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface>
LOCA::TurningPoint::MinimallyAugmented::Constraint::
clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new Constraint(*this, type));
}

int
LOCA::TurningPoint::MinimallyAugmented::Constraint::
numConstraints() const
{
  return 1;
}

void
LOCA::TurningPoint::MinimallyAugmented::Constraint::
setX(const NOX::Abstract::Vector& y)
{
  grpPtr->setX(y);
  isValidConstraints = false;
  isValidDX = false;
}

void
LOCA::TurningPoint::MinimallyAugmented::Constraint::
setParam(int paramID, double val)
{
  grpPtr->setParam(paramID, val);
  isValidConstraints = false;
  isValidDX = false;
}

void
LOCA::TurningPoint::MinimallyAugmented::Constraint::
setParams(const vector<int>& paramIDs, 
	  const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  grpPtr->setParamsMulti(paramIDs, vals);
  isValidConstraints = false;
  isValidDX = false;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::Constraint::
computeConstraints()
{
  if (isValidConstraints)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::TurningPoint::MinimallyAugmented::Constraint::computeConstraints()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // Compute J
  status = grpPtr->computeJacobian();
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							   finalStatus,
							   callingFunction);

  // Set up bordered systems
  borderedSolver->setMatrixBlocksMultiVecConstraint(grpPtr, 
						    a_vector, 
						    b_vector, 
						    Teuchos::null);

  // Create RHS
  NOX::Abstract::MultiVector::DenseMatrix one(1,1);
  one(0,0) = dn;

  // Get linear solver parameters
  Teuchos::RefCountPtr<Teuchos::ParameterList> linear_solver_params =
    parsedParams->getSublist("Linear Solver");

  // Compute sigma_1 and right null vector v
  NOX::Abstract::MultiVector::DenseMatrix s1(1,1);
  status = borderedSolver->initForSolve();
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							   finalStatus,
							   callingFunction);
  status = borderedSolver->applyInverse(*linear_solver_params, 
					NULL, 
					&one, 
					*v_vector, 
					s1);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							   finalStatus,
							   callingFunction);

  // Compute sigma_2 and left null vector w
  NOX::Abstract::MultiVector::DenseMatrix s2(1,1);
  if (!isSymmetric) {
    status = borderedSolver->initForTransposeSolve();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
    status = borderedSolver->applyInverseTranspose(*linear_solver_params, 
						   NULL, 
						   &one, 
						   *w_vector, 
						   s2);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);

  }
  else {
    *w_vector = *v_vector;
    s2.assign(s1);
  }
  
  // Compute sigma = -w^T*J*v
  status = grpPtr->applyJacobianMultiVector(*v_vector, *Jv_vector);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							   finalStatus,
							   callingFunction);
  Jv_vector->multiply(-1.0, *w_vector, constraints);

  // Scale sigma
//   double w_norm = (*w_vector)[0].norm();
//   double v_norm = (*v_vector)[0].norm();
//   sigma_scale = w_norm*v_norm;
  sigma_scale = dn;
  constraints.scale(1.0/sigma_scale);

  if (globalData->locaUtils->isPrintType(NOX::Utils::OuterIteration)) {
    globalData->locaUtils->out() << 
      "\n\tEstimate for singularity of Jacobian (sigma1) = " << 
      globalData->locaUtils->sciformat(s1(0,0));
    globalData->locaUtils->out() << 
      "\n\tEstimate for singularity of Jacobian (sigma2) = " << 
      globalData->locaUtils->sciformat(s2(0,0));
    globalData->locaUtils->out() << 
      "\n\tEstimate for singularity of Jacobian (sigma) = " << 
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

    a_vector->scale(std::sqrt(dn) / (*a_vector)[0].norm());
    b_vector->scale(std::sqrt(dn) / (*b_vector)[0].norm());
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::Constraint::
computeDX()
{
  if (isValidDX)
    return NOX::Abstract::Group::Ok;

  string callingFunction = 
    "LOCA::TurningPoint::MinimallyAugmented::Constraint::computeDX()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // Compute sigma, w and v if necessary
  if (!isValidConstraints) {
    status = computeConstraints();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  // Compute -(w^T*J*v)_x
  status = grpPtr->computeDwtJnDx((*w_vector)[0], (*v_vector)[0], 
				  (*sigma_x)[0]);
  finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  sigma_x->scale(-1.0/sigma_scale);

  isValidDX = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::Constraint::
computeDP(const vector<int>& paramIDs, 
	  NOX::Abstract::MultiVector::DenseMatrix& dgdp, 
	  bool isValidG)
{
  string callingFunction = 
    "LOCA::TurningPoint::MinimallyAugmented::Constraint::computeDP()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // Compute sigma, w and v if necessary
  if (!isValidConstraints) {
    status = computeConstraints();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  // Compute -(w^T*J*v)_p
  status = grpPtr->computeDwtJnDp(paramIDs, (*w_vector)[0], (*v_vector)[0], 
				  dgdp, false);
  finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  dgdp.scale(-1.0/sigma_scale);

  // Set the first column of dgdp
  dgdp(0,0) = constraints(0,0);

  return finalStatus;
}

bool
LOCA::TurningPoint::MinimallyAugmented::Constraint::
isConstraints() const
{
  return isValidConstraints;
}

bool
LOCA::TurningPoint::MinimallyAugmented::Constraint::
isDX() const
{
  return isValidDX;
}

const NOX::Abstract::MultiVector::DenseMatrix&
LOCA::TurningPoint::MinimallyAugmented::Constraint::
getConstraints() const
{
  return constraints;
}

const NOX::Abstract::MultiVector*
LOCA::TurningPoint::MinimallyAugmented::Constraint::
getDX() const
{
  return sigma_x.get();
}

bool
LOCA::TurningPoint::MinimallyAugmented::Constraint::
isDXZero() const
{
  return false;
}

void
LOCA::TurningPoint::MinimallyAugmented::Constraint::
postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  if (stepStatus == LOCA::Abstract::Iterator::Successful && 
      updateVectorsEveryContinuationStep) {
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
      globalData->locaUtils->out() << 
      "\n\tUpdating null vectors for the next continuation step" << std::endl;
    }
    *a_vector = *w_vector;
    *b_vector = *v_vector;

    a_vector->scale(std::sqrt(dn) / (*a_vector)[0].norm());
    b_vector->scale(std::sqrt(dn) / (*b_vector)[0].norm());
  }
}
