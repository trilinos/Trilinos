// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
    TEUCHOS_TEST_FOR_EXCEPTION(ret != Converged || ret != Unconverged,
		       std::invalid_argument,
		       "Invalid ReturnType enum value " << ret << ".  "
		       "Valid values are Converged=" << Converged << " and "
		       "Unconverged=" << Unconverged << ".");
    return ret;
  }

  int 
  InnerSolveResult::requireNonNegInt (const int k) 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(k < 0, std::invalid_argument, 
		       "The given integer argument k=" << k 
		       << " must be nonnegative.");
    return k;
  }

} // namespace Belos
