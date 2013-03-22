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

#include "LOCA_StatusTest_MaxIters.H" // class definition
#include "LOCA_StatusTest_Abstract.H"

//#include "LOCA_Stepper.H"
#include "LOCA_Abstract_Iterator.H"

// FIXME remove these headers?
#include "NOX_Utils.H"
#include "LOCA_GlobalData.H"

LOCA::StatusTest::MaxIters::
MaxIters(int maxIterations,
         bool return_failed_on_max_steps,
         const Teuchos::RCP<const LOCA::GlobalData> globalDataPtr ) :
  maxiters(maxIterations),
  return_failed_on_max_steps_(return_failed_on_max_steps),
  niters(0),
  status(LOCA::StatusTest::Unevaluated)
{
  if ( globalDataPtr.is_valid_ptr() && !globalDataPtr.is_null() )
    globalDataPtr_ = globalDataPtr;

  if (maxiters < 0)
  {
    if ( globalDataPtr_.is_valid_ptr() && !globalDataPtr_.is_null() )
        globalDataPtr_->locaUtils->err() << "LOCA::StatusTest::MaxIters - must choose a number greater than or equal to zero" << std::endl;
    else
        // This will spit out the error message NUMPROC times. -- Without locaUtils, there's nothing we can do..
        std::cerr << "LOCA::StatusTest::MaxIters - must choose a number greater than or equal to zero" << std::endl;
    throw "LOCA Error";
  }
}

LOCA::StatusTest::MaxIters::~MaxIters()
{
}

LOCA::StatusTest::StatusType LOCA::StatusTest::MaxIters::
//checkStatus(const LOCA::Stepper& stepper,
checkStatus(const LOCA::Abstract::Iterator& stepper,
	    LOCA::StatusTest::CheckType checkType)
{
  switch (checkType)
  {
  case LOCA::StatusTest::Complete:
  case LOCA::StatusTest::Minimal:
    niters = stepper.getNumTotalSteps();
    if (niters >= maxiters)
      status = return_failed_on_max_steps_ ? LOCA::StatusTest::Failed : LOCA::StatusTest::Finished;
    else
      status = LOCA::StatusTest::NotFinished;
    break;

  case LOCA::StatusTest::None:
  default:
    niters = -1;
    status = LOCA::StatusTest::Unevaluated;
    break;
  }

  return status;
}

LOCA::StatusTest::StatusType LOCA::StatusTest::MaxIters::
getStatus() const
{
  return status;
}

std::ostream& LOCA::StatusTest::MaxIters::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "Number of Iterations = " << niters << " < " << maxiters;
  stream << std::endl;
 return stream;
}

int LOCA::StatusTest::MaxIters::getMaxIters() const
{
  return maxiters;
}

int LOCA::StatusTest::MaxIters::getNumIters() const
{
  return niters;
}
