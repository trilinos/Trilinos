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
#include "LOCA_Abstract_Iterator.H"    // class definition
#include "NOX_Parameter_List.H"

LOCA::Abstract::Iterator::Iterator()
{
}

LOCA::Abstract::Iterator::Iterator(NOX::Parameter::List& p)
{
  reset(p);
}

LOCA::Abstract::Iterator::Iterator(const LOCA::Abstract::Iterator& it) :
  stepNumber(it.stepNumber),
  numFailedSteps(it.numFailedSteps),
  numTotalSteps(it.numTotalSteps),
  maxSteps(it.maxSteps),
  iteratorStatus(it.iteratorStatus)
{}

LOCA::Abstract::Iterator::~Iterator() {}

bool 
LOCA::Abstract::Iterator::reset(NOX::Parameter::List& p) 
{
  stepNumber = 0;
  numFailedSteps = 0;
  numTotalSteps = 0;
  iteratorStatus = LOCA::Abstract::Iterator::NotFinished;

  maxSteps = p.getParameter("Max Steps",100);
  
  return true;
}

LOCA::Abstract::Iterator::IteratorStatus 
LOCA::Abstract::Iterator::getIteratorStatus() const
{
  return iteratorStatus;
}

int 
LOCA::Abstract::Iterator::getStepNumber() const
{
  return stepNumber;
}

int 
LOCA::Abstract::Iterator::getNumFailedSteps() const
{
  return numFailedSteps;
}

int
LOCA::Abstract::Iterator::getNumTotalSteps() const
{
  return numTotalSteps;
}

LOCA::Abstract::Iterator::IteratorStatus 
LOCA::Abstract::Iterator::run() 
{
  iteratorStatus = start();
  if (iteratorStatus == LOCA::Abstract::Iterator::Finished ||
      iteratorStatus == LOCA::Abstract::Iterator::Failed)
    return iteratorStatus;

  iteratorStatus = iterate();
  if (iteratorStatus == LOCA::Abstract::Iterator::Finished ||
      iteratorStatus == LOCA::Abstract::Iterator::Failed)
    return iteratorStatus;

  iteratorStatus = finish();
 
  return iteratorStatus;
}

LOCA::Abstract::Iterator::IteratorStatus 
LOCA::Abstract::Iterator::iterate()
{
  LOCA::Abstract::Iterator::StepStatus stepStatus = 
    LOCA::Abstract::Iterator::Successful;

  while (iteratorStatus == LOCA::Abstract::Iterator::NotFinished) {

    stepStatus = preprocess(stepStatus);

    stepStatus = compute(stepStatus);

    stepStatus = postprocess(stepStatus);

    ++numTotalSteps;
    if (stepStatus =  LOCA::Abstract::Iterator::Successful)
      ++stepNumber;
    else
      ++numFailedSteps;

    iteratorStatus = stop(stepStatus);
  }

  return iteratorStatus;
}

LOCA::Abstract::Iterator::IteratorStatus 
LOCA::Abstract::Iterator::stop(StepStatus stepStatus)
{
  if (numTotalSteps >= maxSteps)
    return LOCA::Abstract::Iterator::Finished;
  else
    return LOCA::Abstract::Iterator::NotFinished;
}
