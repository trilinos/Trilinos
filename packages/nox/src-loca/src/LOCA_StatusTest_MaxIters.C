// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
