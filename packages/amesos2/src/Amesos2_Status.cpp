/**
  \file   Amesos2_Status.cpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Sat Jan 16 09:44:07 2010

  \brief  Implementation for Amesos2::Status
*/
#ifndef AMESOS2_STATUS_CPP
#define AMESOS2_STATUS_CPP

#include "Amesos2_Status.hpp"

namespace Amesos {


void Status::setStatusParameters(
  const Teuchos::RCP<Teuchos::ParameterList> & parameterList) {

  // some verbose output:
  // 0 - no output at all
  // 1 - output as specified by other parameters
  // 2 - all possible output
  if( parameterList->isParameter("OutputLevel") )
    verbose_ = parameterList->get<int>("OutputLevel");

  // level of debug output:
  // 0 - no output at all
  // 1 - some debug output - set by some tests upon a test failure
  // >1 - more debug output (unused at this point)
  if( parameterList->isParameter("DebugLevel") )
    debug_ = parameterList->get<int>("DebugLevel");

  // print some timing information (on process 0)
  if( parameterList->isParameter("PrintTiming") )
    printTiming_ = parameterList->get<bool>("PrintTiming");

  // print some statistics (on process 0). Do not include timing
  if( parameterList->isParameter("PrintStatus") )
    printStatus_ = parameterList->get<bool>("PrintStatus");

  // compute norms of some vectors
  if( parameterList->isParameter("ComputeVectorNorms") )
    computeVectorNorms_ = parameterList->get<bool>("ComputeVectorNorms");

  // compute the true residual Ax-b after solution
  if( parameterList->isParameter("ComputeTrueResidual") )
    computeTrueResidual_ = parameterList->get<bool>("ComputeTrueResidual");

}


} // end namespace Amesos

#endif	// AMESOS2_STATUS_CPP
