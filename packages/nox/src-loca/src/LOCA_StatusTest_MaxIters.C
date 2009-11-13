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

#include "LOCA_StatusTest_MaxIters.H" // class definition
#include "LOCA_StatusTest_Abstract.H"

#include "LOCA_Stepper.H"

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

  if (maxiters < 1) 
  {
    globalDataPtr_->locaUtils->err() << "LOCA::StatusTest::MaxIters - must choose a number greater than zero" << endl;
    throw "LOCA Error";
  }
}

LOCA::StatusTest::MaxIters::~MaxIters()
{
}

LOCA::StatusTest::StatusType LOCA::StatusTest::MaxIters::
checkStatus(const LOCA::Stepper& stepper,
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

ostream& LOCA::StatusTest::MaxIters::print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "Number of Iterations = " << niters << " < " << maxiters;
  stream << endl;
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
