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
LOCA::Abstract::Iterator::stop(LOCA::Abstract::Iterator::StepStatus /* stepStatus */)
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
