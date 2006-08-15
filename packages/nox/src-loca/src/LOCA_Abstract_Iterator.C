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
#include "LOCA_Abstract_Iterator.H"    // class definition
#include "Teuchos_ParameterList.hpp"

LOCA::Abstract::Iterator::Iterator() :
  stepNumber(0),
  numFailedSteps(0),
  numTotalSteps(0),
  maxSteps(100),
  iteratorStatus(LOCA::Abstract::Iterator::NotFinished)
{
}

LOCA::Abstract::Iterator::Iterator(Teuchos::ParameterList& p) :
  stepNumber(0),
  numFailedSteps(0),
  numTotalSteps(0),
  maxSteps(100),
  iteratorStatus(LOCA::Abstract::Iterator::NotFinished)
{
  resetIterator(p);
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
LOCA::Abstract::Iterator::resetIterator(Teuchos::ParameterList& p) 
{
  stepNumber = 0;
  numFailedSteps = 0;
  numTotalSteps = 0;
  iteratorStatus = LOCA::Abstract::Iterator::NotFinished;

  maxSteps = p.get("Max Steps",100);
  
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
  if (iteratorStatus == LOCA::Abstract::Iterator::Failed)
    return iteratorStatus;

  stepNumber++;

  iteratorStatus = iterate();

  iteratorStatus = finish(iteratorStatus);
 
  return iteratorStatus;
}

LOCA::Abstract::Iterator::IteratorStatus 
LOCA::Abstract::Iterator::iterate()
{
  LOCA::Abstract::Iterator::StepStatus stepStatus = 
    LOCA::Abstract::Iterator::Successful;
  LOCA::Abstract::Iterator::StepStatus preStatus;
  LOCA::Abstract::Iterator::StepStatus compStatus;
  LOCA::Abstract::Iterator::StepStatus postStatus;

  iteratorStatus = stop(stepStatus);

  while (iteratorStatus == LOCA::Abstract::Iterator::NotFinished) {

    preStatus = preprocess(stepStatus);

    compStatus = compute(preStatus);

    postStatus = postprocess(compStatus);

    stepStatus = computeStepStatus(preStatus, compStatus, postStatus);

    ++numTotalSteps;
    if (stepStatus ==  LOCA::Abstract::Iterator::Successful)
      ++stepNumber;
    else
      ++numFailedSteps;

    if (iteratorStatus != LOCA::Abstract::Iterator::Failed)
      iteratorStatus = stop(stepStatus);
  }

  return iteratorStatus;
}

LOCA::Abstract::Iterator::IteratorStatus 
LOCA::Abstract::Iterator::stop(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  if (numTotalSteps >= maxSteps)
    return LOCA::Abstract::Iterator::Finished;
  else
    return LOCA::Abstract::Iterator::NotFinished;
}

void
LOCA::Abstract::Iterator::setLastIteration() 
{
  iteratorStatus = LOCA::Abstract::Iterator::LastIteration;
  return;
}

bool
LOCA::Abstract::Iterator::isLastIteration() 
{
  return (iteratorStatus == LOCA::Abstract::Iterator::LastIteration);
}

LOCA::Abstract::Iterator::StepStatus 
LOCA::Abstract::Iterator::computeStepStatus(
		       LOCA::Abstract::Iterator::StepStatus preStatus, 
		       LOCA::Abstract::Iterator::StepStatus compStatus,
		       LOCA::Abstract::Iterator::StepStatus postStatus)
{
  bool haveProvisional = false;
  bool haveUnsuccessful = false;

  if (preStatus == LOCA::Abstract::Iterator::Provisional ||
      compStatus == LOCA::Abstract::Iterator::Provisional ||
      postStatus == LOCA::Abstract::Iterator::Provisional) {
    haveProvisional = true;
  }

  if (preStatus == LOCA::Abstract::Iterator::Unsuccessful ||
      compStatus == LOCA::Abstract::Iterator::Unsuccessful ||
      postStatus == LOCA::Abstract::Iterator::Unsuccessful) {
    haveUnsuccessful = true;
  }

  if (haveProvisional && haveUnsuccessful) {
    iteratorStatus = LOCA::Abstract::Iterator::Failed;
    return LOCA::Abstract::Iterator::Unsuccessful;
  }
  else if (haveUnsuccessful)
    return LOCA::Abstract::Iterator::Unsuccessful;
  else 
    return LOCA::Abstract::Iterator::Successful;
}
