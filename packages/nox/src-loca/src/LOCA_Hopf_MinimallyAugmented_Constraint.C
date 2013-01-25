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

#include "LOCA_Hopf_MinimallyAugmented_Constraint.H"
#include "LOCA_Hopf_MinimallyAugmented_AbstractGroup.H"
#include "LOCA_BorderedSolver_AbstractStrategy.H"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Factory.H"
#include "NOX_Utils.H"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_Hopf_ComplexMultiVector.H"
#include "LOCA_BorderedSolver_ComplexOperator.H"

LOCA::Hopf::MinimallyAugmented::Constraint::
Constraint(
    const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
    const Teuchos::RCP<Teuchos::ParameterList>& hpfParams,
    const Teuchos::RCP<LOCA::Hopf::MinimallyAugmented::AbstractGroup>& g,
    bool is_symmetric,
    const NOX::Abstract::Vector& a_real,
    const NOX::Abstract::Vector& a_imag,
    const NOX::Abstract::Vector* b_real,
    const NOX::Abstract::Vector* b_imag,
    int bif_param,
    double freq) :
  globalData(global_data),
  parsedParams(topParams),
  hopfParams(hpfParams),
  grpPtr(g),
  a_vector(a_real.createMultiVector(2, NOX::ShapeCopy)),
  b_vector(),
  w_vector(a_real.createMultiVector(2, NOX::ShapeCopy)),
  v_vector(a_real.createMultiVector(2, NOX::ShapeCopy)),
  Cv_vector(a_real.createMultiVector(2, NOX::ShapeCopy)),
  sigma_x(a_real.createMultiVector(2, NOX::ShapeCopy)),
  constraints(2, 1),
  borderedSolver(),
  dn(static_cast<double>(a_vector->length())),
  sigma_scale(1.0),
  isSymmetric(is_symmetric),
  isValidConstraints(false),
  isValidDX(false),
  bifParamID(1),
  omega(freq),
  updateVectorsEveryContinuationStep(true),
  updateVectorsEveryIteration(false)
{
  // Instantiate bordered solvers
  borderedSolver = 
    globalData->locaFactory->createBorderedSolverStrategy(parsedParams,
							  hopfParams);

  (*a_vector)[0] = a_real;
  (*a_vector)[1] = a_imag;
  if (!isSymmetric) {
    b_vector = b_real->createMultiVector(2, NOX::ShapeCopy);
    (*b_vector)[0] = *b_real;
    (*b_vector)[1] = *b_imag;
  }
  else {
    b_vector = a_vector->clone(NOX::DeepCopy);
  }

  // Options
  updateVectorsEveryContinuationStep = 
    hopfParams->get("Update Null Vectors Every Continuation Step", 
		    true);
  updateVectorsEveryIteration = 
    hopfParams->get("Update Null Vectors Every Nonlinear Iteration", 
		    false);
}

LOCA::Hopf::MinimallyAugmented::Constraint::
Constraint(const LOCA::Hopf::MinimallyAugmented::Constraint& source, 
	   NOX::CopyType type) : 
  globalData(source.globalData),
  parsedParams(source.parsedParams),
  hopfParams(source.hopfParams),
  grpPtr(Teuchos::null),
  a_vector(source.a_vector->clone(type)),
  b_vector(source.b_vector->clone(type)),
  w_vector(source.w_vector->clone(type)),
  v_vector(source.v_vector->clone(type)),
  Cv_vector(source.Cv_vector->clone(type)),
  sigma_x(source.sigma_x->clone(type)),
  constraints(source.constraints),
  borderedSolver(),
  dn(source.dn),
  sigma_scale(source.sigma_scale),
  isSymmetric(source.isSymmetric),
  isValidConstraints(false),
  isValidDX(false),
  bifParamID(source.bifParamID),
  omega(source.omega),
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
							  hopfParams);

  // We don't explicitly copy the group because the constrained group
  // will do that
}

LOCA::Hopf::MinimallyAugmented::Constraint::
~Constraint()
{
}

void
LOCA::Hopf::MinimallyAugmented::Constraint::
setGroup(const Teuchos::RCP<LOCA::Hopf::MinimallyAugmented::AbstractGroup>& g)
{
  grpPtr = g;
}

void
LOCA::Hopf::MinimallyAugmented::Constraint::
setFrequency(double freq)
{
  omega = freq;
  isValidConstraints = false;
  isValidDX = false;
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::Hopf::MinimallyAugmented::Constraint::
getLeftNullVecReal() const
{
  return Teuchos::rcp(&(*w_vector)[0], false);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::Hopf::MinimallyAugmented::Constraint::
getLeftNullVecImag() const
{
  return Teuchos::rcp(&(*w_vector)[1], false);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::Hopf::MinimallyAugmented::Constraint::
getRightNullVecReal() const
{
  return Teuchos::rcp(&(*v_vector)[0], false);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::Hopf::MinimallyAugmented::Constraint::
getRightNullVecImag() const
{
  return Teuchos::rcp(&(*v_vector)[1], false);
}

double
LOCA::Hopf::MinimallyAugmented::Constraint::
getSigmaReal() const
{
  return constraints(0,0);
}

double
LOCA::Hopf::MinimallyAugmented::Constraint::
getSigmaImag() const
{
  return constraints(1,0);
}

void
LOCA::Hopf::MinimallyAugmented::Constraint::
copy(const LOCA::MultiContinuation::ConstraintInterface& src)
{
  const LOCA::Hopf::MinimallyAugmented::Constraint& source = 
  dynamic_cast<const LOCA::Hopf::MinimallyAugmented::Constraint&>(src);

  if (this != &source) {
    globalData = source.globalData;
    parsedParams = source.parsedParams;
    hopfParams = source.hopfParams;
    *a_vector = *(source.a_vector);
    *b_vector = *(source.b_vector);
    *w_vector = *(source.w_vector);
    *v_vector = *(source.v_vector);
    *Cv_vector = *(source.Cv_vector);
    *sigma_x = *(source.sigma_x);
    constraints.assign(source.constraints);
    dn = source.dn;
    sigma_scale = source.sigma_scale;
    isSymmetric = source.isSymmetric;
    isValidConstraints = source.isValidConstraints;
    isValidDX = source.isValidDX;
    bifParamID = source.bifParamID;
    omega = source.omega;
    updateVectorsEveryContinuationStep = 
      source.updateVectorsEveryContinuationStep;
    updateVectorsEveryIteration = 
      source.updateVectorsEveryIteration;

    // Instantiate bordered solvers
    borderedSolver = 
      globalData->locaFactory->createBorderedSolverStrategy(parsedParams,
							    hopfParams);

    // We don't explicitly copy the group because the constrained group
    // will do that
  }
}

Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
LOCA::Hopf::MinimallyAugmented::Constraint::
clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new Constraint(*this, type));
}

int
LOCA::Hopf::MinimallyAugmented::Constraint::
numConstraints() const
{
  return 2;
}

void
LOCA::Hopf::MinimallyAugmented::Constraint::
setX(const NOX::Abstract::Vector& y)
{
  grpPtr->setX(y);
  isValidConstraints = false;
  isValidDX = false;
}

void
LOCA::Hopf::MinimallyAugmented::Constraint::
setParam(int paramID, double val)
{
  grpPtr->setParam(paramID, val);
  isValidConstraints = false;
  isValidDX = false;
}

void
LOCA::Hopf::MinimallyAugmented::Constraint::
setParams(const std::vector<int>& paramIDs, 
	  const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  grpPtr->setParamsMulti(paramIDs, vals);
  isValidConstraints = false;
  isValidDX = false;
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::Constraint::
computeConstraints()
{
  if (isValidConstraints)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction = 
    "LOCA::Hopf::MinimallyAugmented::Constraint::computeConstraints()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  // Compute C
  status = grpPtr->computeComplex(omega);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							   finalStatus,
							   callingFunction);

  // Compute A and B blocks
  Teuchos::RCP<LOCA::Hopf::ComplexMultiVector> A = 
    Teuchos::rcp(new LOCA::Hopf::ComplexMultiVector(globalData, 
						    (*a_vector)[0], 2));
  (*(A->getRealMultiVec()))[0] = (*a_vector)[0];
  (*(A->getImagMultiVec()))[0] = (*a_vector)[1];
  (*(A->getRealMultiVec()))[1] = (*a_vector)[1];
  (*(A->getImagMultiVec()))[1] = (*a_vector)[0];
  (*(A->getRealMultiVec()))[1].scale(-1.0);

  Teuchos::RCP<LOCA::Hopf::ComplexMultiVector> B = 
    Teuchos::rcp(new LOCA::Hopf::ComplexMultiVector(globalData, 
						    (*b_vector)[0], 2));
  (*(B->getRealMultiVec()))[0] = (*b_vector)[0];
  (*(B->getImagMultiVec()))[0] = (*b_vector)[1];
  (*(B->getRealMultiVec()))[1] = (*b_vector)[1];
  (*(B->getImagMultiVec()))[1] = (*b_vector)[0];
  (*(B->getRealMultiVec()))[1].scale(-1.0);

  // Set up bordered systems
  Teuchos::RCP<const LOCA::BorderedSolver::ComplexOperator> op =
    Teuchos::rcp(new  LOCA::BorderedSolver::ComplexOperator(grpPtr, omega));
  borderedSolver->setMatrixBlocksMultiVecConstraint(op, A, B, 
						    Teuchos::null);

  // Create RHS
  NOX::Abstract::MultiVector::DenseMatrix one(2,1);
  one(0,0) = dn;
  one(1,0) = 0.0;

  // Get linear solver parameters
  Teuchos::RCP<Teuchos::ParameterList> linear_solver_params =
    parsedParams->getSublist("Linear Solver");

  // Compute sigma_1 and right null vector v
  NOX::Abstract::MultiVector::DenseMatrix s1(2,1);
  Teuchos::RCP<LOCA::Hopf::ComplexMultiVector> V = 
    Teuchos::rcp(new LOCA::Hopf::ComplexMultiVector(globalData, 
						    (*v_vector)[0], 1));
  (*(V->getRealMultiVec()))[0] = (*v_vector)[0];
  (*(V->getImagMultiVec()))[0] = (*v_vector)[1];
  status = borderedSolver->initForSolve();
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							   finalStatus,
							   callingFunction);
  status = borderedSolver->applyInverse(*linear_solver_params, 
					NULL, 
					&one, 
					*V, 
					s1);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							   finalStatus,
							   callingFunction);
  (*v_vector)[0] = (*(V->getRealMultiVec()))[0];
  (*v_vector)[1] = (*(V->getImagMultiVec()))[0];

//   Teuchos::RCP<NOX::Abstract::MultiVector> rV = 
//     v_vector->clone(NOX::ShapeCopy);
//   NOX::Abstract::MultiVector::DenseMatrix rs1(2,1);
//   rV->init(0.0);
//   rs1.putScalar(0.0);
  
//   grpPtr->applyComplex((*v_vector)[0], (*v_vector)[1], (*rV)[0], (*rV)[1]);
//   (*rV)[0].update(s1(0,0), (*a_vector)[0], -s1(1,0), (*a_vector)[1], 1.0);
//   (*rV)[1].update(s1(0,0), (*a_vector)[1],  s1(1,0), (*a_vector)[0], 1.0);

//   rs1(0,0) = (*b_vector)[0].innerProduct((*v_vector)[0]) + 
//     (*b_vector)[1].innerProduct((*v_vector)[1]);
//   rs1(1,0) = (*b_vector)[0].innerProduct((*v_vector)[1]) - 
//     (*b_vector)[1].innerProduct((*v_vector)[0]);
//   rs1 -= one;

//   std::cout << "rV = " << std::endl;
//   rV->print(cout);

//   std::cout << "rs1 = " << std::endl;
//   rs1.print(cout);

//   std::cout << "checking error..." << std::endl;
//   Teuchos::RCP<NOX::Abstract::MultiVector> rV = 
//     V->clone(NOX::ShapeCopy);
//   NOX::Abstract::MultiVector::DenseMatrix rs1(2,1);
//   borderedSolver->apply(*V, s1, *rV, rs1);
//   rs1 -= one;
//   std::cout << "rV->norm() = " << (*rV)[0].norm() << std::endl;
//   std::cout << "rs1.norm() = " << rs1.normInf() << std::endl;

  // Compute sigma_2 and left null vector w
  NOX::Abstract::MultiVector::DenseMatrix s2(2,1);
  Teuchos::RCP<LOCA::Hopf::ComplexMultiVector> W = 
    Teuchos::rcp(new LOCA::Hopf::ComplexMultiVector(globalData, 
						    (*w_vector)[0], 1));
  (*(W->getRealMultiVec()))[0] = (*w_vector)[0];
  (*(W->getImagMultiVec()))[0] = (*w_vector)[1];
  if (!isSymmetric) {
    status = borderedSolver->initForTransposeSolve();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
    status = borderedSolver->applyInverseTranspose(*linear_solver_params, 
						   NULL, 
						   &one, 
						   *W, 
						   s2);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);

    (*w_vector)[0] = (*(W->getRealMultiVec()))[0];
    (*w_vector)[1] = (*(W->getImagMultiVec()))[0];

  }
  else {
    *w_vector = *v_vector;
    s2.assign(s1);
  }

//   Teuchos::RCP<NOX::Abstract::MultiVector> rW = 
//     v_vector->clone(NOX::ShapeCopy);
//   NOX::Abstract::MultiVector::DenseMatrix rs2(2,1);
//   rW->init(0.0);
//   rs2.putScalar(0.0);
  
//   grpPtr->applyComplexTranspose((*w_vector)[0], (*w_vector)[1], 
// 				(*rW)[0], (*rW)[1]);
//   (*rW)[0].update(s2(0,0), (*b_vector)[0], -s2(1,0), (*b_vector)[1], 1.0);
//   (*rW)[1].update(s2(0,0), (*b_vector)[1],  s2(1,0), (*b_vector)[0], 1.0);

//   rs2(0,0) = (*a_vector)[0].innerProduct((*w_vector)[0]) + 
//     (*a_vector)[1].innerProduct((*w_vector)[1]);
//   rs2(1,0) = (*a_vector)[0].innerProduct((*w_vector)[1]) - 
//     (*a_vector)[1].innerProduct((*w_vector)[0]);
//   rs2 -= one;

//   std::cout << "rW = " << std::endl;
//   rW->print(cout);

//   std::cout << "rs2 = " << std::endl;
//   rs2.print(cout);

//   Teuchos::RCP<NOX::Abstract::MultiVector> rW = 
//     W->clone(NOX::ShapeCopy);
//   NOX::Abstract::MultiVector::DenseMatrix rs2(2,1);
//   borderedSolver->applyTranspose(*W, s2, *rW, rs2);
//   rs2 -= one;
//   std::cout << "rW->norm() = " << (*rW)[0].norm() << std::endl;
//   std::cout << "rs2.norm() = " << rs2.normInf() << std::endl;
  
  // Compute sigma = -w^T*J*v
  status = grpPtr->applyComplex((*v_vector)[0], (*v_vector)[1],
				(*Cv_vector)[0], (*Cv_vector)[1]);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							   finalStatus,
							   callingFunction);
  NOX::Abstract::MultiVector::DenseMatrix tmp(2,2);
  Cv_vector->multiply(-1.0, *w_vector, tmp);
  constraints(0,0) = tmp(0,0) + tmp(1,1);
  constraints(1,0) = tmp(0,1) - tmp(1,0);

  // Scale sigma
  sigma_scale = dn;
  constraints.scale(1.0/sigma_scale);

  if (globalData->locaUtils->isPrintType(NOX::Utils::OuterIteration)) {
    globalData->locaUtils->out() << 
      "\n\tEstimate for singularity of Complex Jacobian (sigma1) =\n\t\t" << 
      globalData->locaUtils->sciformat(s1(0,0));
    if (s1(1,0) > 0.0) 
      globalData->locaUtils->out() << " + i ";
    else
      globalData->locaUtils->out() << " - i ";
    globalData->locaUtils->out() <<
      globalData->locaUtils->sciformat(std::fabs(s1(1,0)));
    globalData->locaUtils->out() << 
      "\n\tEstimate for singularity of Complex Jacobian (sigma2) =\n\t\t" << 
      globalData->locaUtils->sciformat(s2(0,0));
    if (s2(1,0) > 0.0) 
      globalData->locaUtils->out() << " + i ";
    else
      globalData->locaUtils->out() << " - i ";
    globalData->locaUtils->out() << 
      globalData->locaUtils->sciformat(std::fabs(s2(1,0)));
    globalData->locaUtils->out() << 
      "\n\tEstimate for singularity of Complex Jacobian (sigma ) =\n\t\t" << 
      globalData->locaUtils->sciformat(constraints(0,0));
    if (constraints(1,0) > 0.0)
      globalData->locaUtils->out() << " + i ";
    else
      globalData->locaUtils->out() << " - i ";
    globalData->locaUtils->out() <<
      globalData->locaUtils->sciformat(std::fabs(constraints(1,0))) << std::endl;
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

    double a1n = (*a_vector)[0].norm();
    double a2n = (*a_vector)[1].norm();
    double b1n = (*b_vector)[0].norm();
    double b2n = (*b_vector)[1].norm();
    a_vector->scale(std::sqrt(dn) / std::sqrt(a1n*a1n + a2n*a2n));
    b_vector->scale(std::sqrt(dn) / std::sqrt(b1n*b1n + b2n*b2n));
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::Constraint::
computeDX()
{
  if (isValidDX)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction = 
    "LOCA::Hopf::MinimallyAugmented::Constraint::computeDX()";
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
  status = grpPtr->computeDwtCeDx((*w_vector)[0], (*w_vector)[1], 
				  (*v_vector)[0], (*v_vector)[1],
				  omega,
				  (*sigma_x)[0], (*sigma_x)[1]);
  finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  sigma_x->scale(-1.0/sigma_scale);

  isValidDX = true;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::Constraint::
computeDP(const std::vector<int>& paramIDs, 
	  NOX::Abstract::MultiVector::DenseMatrix& dgdp, 
	  bool isValidG)
{
  std::string callingFunction = 
    "LOCA::Hopf::MinimallyAugmented::Constraint::computeDP()";
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
  NOX::Abstract::MultiVector::DenseMatrix dgdp_real(Teuchos::View, dgdp,
						    1, paramIDs.size()+1,
						    0, 0);
  NOX::Abstract::MultiVector::DenseMatrix dgdp_imag(Teuchos::View, dgdp,
						    1, paramIDs.size()+1,
						    1, 0);
  status = grpPtr->computeDwtCeDp(paramIDs, 
				  (*w_vector)[0], (*w_vector)[1],
				  (*v_vector)[0], (*v_vector)[1],
				  omega, 
				  dgdp_real, dgdp_imag, false);
  finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  dgdp.scale(-1.0/sigma_scale);

  // Set the first column of dgdp
  dgdp(0,0) = constraints(0,0);
  dgdp(1,0) = constraints(1,0);

  return finalStatus;
}

bool
LOCA::Hopf::MinimallyAugmented::Constraint::
isConstraints() const
{
  return isValidConstraints;
}

bool
LOCA::Hopf::MinimallyAugmented::Constraint::
isDX() const
{
  return isValidDX;
}

const NOX::Abstract::MultiVector::DenseMatrix&
LOCA::Hopf::MinimallyAugmented::Constraint::
getConstraints() const
{
  return constraints;
}

const NOX::Abstract::MultiVector*
LOCA::Hopf::MinimallyAugmented::Constraint::
getDX() const
{
  return sigma_x.get();
}

bool
LOCA::Hopf::MinimallyAugmented::Constraint::
isDXZero() const
{
  return false;
}

void
LOCA::Hopf::MinimallyAugmented::Constraint::
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

    double a1n = (*a_vector)[0].norm();
    double a2n = (*a_vector)[1].norm();
    double b1n = (*b_vector)[0].norm();
    double b2n = (*b_vector)[1].norm();
    a_vector->scale(std::sqrt(dn) / std::sqrt(a1n*a1n + a2n*a2n));
    b_vector->scale(std::sqrt(dn) / std::sqrt(b1n*b1n + b2n*b2n));
  }
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::Constraint::
computeDOmega(NOX::Abstract::MultiVector::DenseMatrix& domega)
{
  std::string callingFunction = 
    "LOCA::Hopf::MinimallyAugmented::Constraint::computeDOmega()";
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

  // Compute mass matrix M
  status = grpPtr->computeShiftedMatrix(0.0, 1.0);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							   finalStatus,
							   callingFunction);

  // Compute M*v
  Teuchos::RCP<NOX::Abstract::MultiVector> tmp_vector =
    v_vector->clone(NOX::ShapeCopy);
  status = grpPtr->applyShiftedMatrixMultiVector(*v_vector, *tmp_vector);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							   finalStatus,
							   callingFunction);

  // Compute u^T*M*v
  NOX::Abstract::MultiVector::DenseMatrix tmp_mat(2,2);
  tmp_vector->multiply(1.0, *w_vector, tmp_mat);

  // Compute domega
  domega(0,0) =   tmp_mat(0,1) - tmp_mat(1,0);
  domega(1,0) = -(tmp_mat(0,0) + tmp_mat(1,1));

  domega.scale(1.0/sigma_scale);

  return finalStatus;
}
