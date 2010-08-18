// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package 
//                  Copyright 2010 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

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
  if( parameterList->isParameter("OutputLevel") ){
    verbose_ = parameterList->get<int>("OutputLevel");
  }

  // level of debug output:
  // 0 - no output at all
  // 1 - some debug output - set by some tests upon a test failure
  // >1 - more debug output (unused at this point)
  if( parameterList->isParameter("DebugLevel") ){
    debug_ = parameterList->get<int>("DebugLevel");
  }

  // print some timing information (on process 0)
  if( parameterList->isParameter("PrintTiming") ){
    printTiming_ = parameterList->get<bool>("PrintTiming");
  }

  // print some statistics (on process 0). Do not include timing
  if( parameterList->isParameter("PrintStatus") ){
    printStatus_ = parameterList->get<bool>("PrintStatus");
  }

  // compute norms of some vectors
  if( parameterList->isParameter("ComputeVectorNorms") ){
    computeVectorNorms_ = parameterList->get<bool>("ComputeVectorNorms");
  }

  // compute the true residual Ax-b after solution
  if( parameterList->isParameter("ComputeTrueResidual") ){
    computeTrueResidual_ = parameterList->get<bool>("ComputeTrueResidual");
  }

}


} // end namespace Amesos

#endif	// AMESOS2_STATUS_CPP
