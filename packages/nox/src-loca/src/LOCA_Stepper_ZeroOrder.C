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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "LOCA_Stepper_ZeroOrder.H"	// class definition

#include "LOCA_Utils.H"                 // Printing utilities
#include "LOCA_Solver_Generic.H"        // class member

using namespace LOCA;
using namespace LOCA::Stepper;

ZeroOrder::ZeroOrder(Solver::Generic& s) :
  Generic(s)
{
  reset(s);
}

ZeroOrder::~ZeroOrder() 
{

}

bool ZeroOrder::reset(Solver::Generic& s) 
{
  resetGenericMembers(s);
  printInitializationInfo();
  return true;
}

bool ZeroOrder::resetFromFailedStep()
{
  // Not implemented yet!
  return false;
}

StatusType ZeroOrder::getStatus()
{
  return status;
}

StatusType ZeroOrder::solve()
{
  while (status == Unconverged) {
    status = step();
  }
  return status;
}

StatusType ZeroOrder::step()
{
  status = NOX::StatusTest::Unconverged;

  solverPtr->reset(LOCA::Solver::PreviousSolution);

  if (stepNumber != 0) {
    prevStepSize = curStepSize;
    prevValue = curValue;
    curStepSize = computeStepSize(solverStatus);
    curValue += curStepSize;
  }      
  
  printStartStep();
  
  conParams.setValue(conParamID, curValue);
  solverPtr->setParams(conParams);
  solverStatus = Unconverged;
  solverStatus = solverPtr->solve();
  
  printEndStep(solverStatus);
  
  if (solverStatus == Failed) {
    numFailedSteps += 1;
    resetFromFailedStep();
  }
  else  {
    stepNumber += 1;
  }
  
  numTotalSteps += 1;

  status = checkStepperStatus();

  return status;
}

const Abstract::Group& ZeroOrder::getSolutionGroup() const
{
  return solverPtr->getSolutionGroup();
}

const Abstract::Group& ZeroOrder::getPreviousSolutionGroup() const
{
  // We should save the soln at the previous time step.
  // The following is wrong - it is the prev nonlinear step.
  return solverPtr->getPreviousSolutionGroup();
}

const NOX::Parameter::List& ZeroOrder::getParameterList() const
{
  return solverPtr->getParameterList();
}
