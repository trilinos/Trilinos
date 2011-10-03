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

ostream& NOX::StatusTest::NStep::print(ostream& stream, int indent) const
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

