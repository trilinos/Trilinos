//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
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
// ************************************************************************
//@HEADER

#include <BelosInnerSolveResult.hpp>

namespace Belos {

  InnerSolveResult::
  InnerSolveResult (const ReturnType theResult,
		    const int theNumRestartCycles,
		    const int theTotalNumIters,
		    const std::map<std::string, double>& theExtraData) :
    result_ (validatedReturnType (theResult)),
    numRestartCycles_ (requireNonNegInt (theNumRestartCycles)),
    totalNumIters_ (requireNonNegInt (theTotalNumIters)),
    extraData_ (theExtraData)
  {}

  InnerSolveResult::
  InnerSolveResult (const ReturnType theResult,
		    const int theNumRestartCycles,
		    const int theTotalNumIters) :
    result_ (validatedReturnType (theResult)),
    numRestartCycles_ (requireNonNegInt (theNumRestartCycles)),
    totalNumIters_ (requireNonNegInt (theTotalNumIters))
  {}

  ReturnType 
  InnerSolveResult::validatedReturnType (const ReturnType ret)
  {
    TEST_FOR_EXCEPTION(ret != Converged || ret != Unconverged,
		       std::invalid_argument,
		       "Invalid ReturnType enum value " << ret << ".  "
		       "Valid values are Converged=" << Converged << " and "
		       "Unconverged=" << Unconverged << ".");
    return ret;
  }

  int 
  InnerSolveResult::requireNonNegInt (const int k) 
  {
    TEST_FOR_EXCEPTION(k < 0, std::invalid_argument, 
		       "The given integer argument k=" << k 
		       << " must be nonnegative.");
    return k;
  }

} // namespace Belos
