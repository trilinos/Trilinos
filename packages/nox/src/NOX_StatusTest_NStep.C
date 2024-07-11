// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Common.H"
#include "NOX_StatusTest_NStep.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"
#include "Teuchos_Assert.hpp"

NOX::StatusTest::NStep::
NStep(int numberOfStepsForConvergence,
      int numberOfNonlinearSolvesInRampingPhase,
      int rampingPhaseNumberOfStepsForConvergence,
      const NOX::Utils* u) :
  status_(Unconverged),
  numberOfStepsForConvergence_(numberOfStepsForConvergence),
  numberOfNonlinearSolvesInRampingPhase_(numberOfNonlinearSolvesInRampingPhase),
  rampingPhaseNumberOfStepsForConvergence_(rampingPhaseNumberOfStepsForConvergence),
  currentNumberOfSteps_(0),
  currentNumberOfNonlinearSolves_(0),
  inRampingPhase_(true)
{
  if (u != NULL)
    utils_ = *u;
}

NOX::StatusTest::StatusType NOX::StatusTest::NStep::
checkStatus(const NOX::Solver::Generic& problem,
        NOX::StatusTest::CheckType checkType)
{

  if (problem.getNumIterations() == 0)
    currentNumberOfNonlinearSolves_ += 1;

  if (currentNumberOfNonlinearSolves_ <= numberOfNonlinearSolvesInRampingPhase_)
    inRampingPhase_ = true;
  else
    inRampingPhase_ = false;

  if (checkType == NOX::StatusTest::None)
  {
    status_ = Unevaluated;
  }
  else
  {
    currentNumberOfSteps_ = problem.getNumIterations();

    if (inRampingPhase_)
      status_ = (currentNumberOfSteps_ >= rampingPhaseNumberOfStepsForConvergence_) ? Converged : Unconverged;
    else
      status_ = (currentNumberOfSteps_ >= numberOfStepsForConvergence_) ? Converged : Unconverged;
  }

  return status_;
}

NOX::StatusTest::StatusType NOX::StatusTest::NStep::getStatus() const
{
  return status_;
}

std::ostream& NOX::StatusTest::NStep::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status_;

  if (inRampingPhase_) {
    stream << "N-Step (in ramping phase " << currentNumberOfNonlinearSolves_
       << "/" << numberOfNonlinearSolvesInRampingPhase_ << "): " << currentNumberOfSteps_;
    if (status_ == Converged)
      stream << " = " << rampingPhaseNumberOfStepsForConvergence_ << std::endl;
    else
      stream << " < " << rampingPhaseNumberOfStepsForConvergence_ << std::endl;
  }
  else {
    stream << "N-Step: " << currentNumberOfSteps_;
    if (status_ == Converged)
      stream << " = " << numberOfStepsForConvergence_ << std::endl;
    else
      stream << " < " << numberOfStepsForConvergence_ << std::endl;
  }

  return stream;
}

