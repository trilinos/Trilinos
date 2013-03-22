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

