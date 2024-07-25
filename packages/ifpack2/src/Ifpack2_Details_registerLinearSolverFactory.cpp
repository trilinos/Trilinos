// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Ifpack2_Details_registerLinearSolverFactory.hpp"
#include "Ifpack2_Details_LinearSolverFactory.hpp"
#include "Ifpack2_ETIHelperMacros.h"

// Define the Ifpack2 ETI macros and Tpetra typedefs that go along
// with them.  This works whether or not ETI is ON, because Bug 6380
// was fixed.
IFPACK2_ETI_MANGLING_TYPEDEFS()

// Macro that registers Ifpack2's LinearSolverFactory for the given
// four template parameters (Scalar = SC, LocalOrdinal = LO,
// GlobalOrdinal = GO, Node = NT).  The macro is local to this file.
//
// NOTE: This macro does NOT do explicit instantiation!  That's why I
// call it LCL_CALL and not LCL_INST.  We are just using the fix for
// Bug 6380 to invoke this class method over the set of enabled
// template parameters.
#define LCL_CALL( SC, LO, GO, NT ) \
  ::Ifpack2::Details::LinearSolverFactory<SC, LO, GO, NT>::registerLinearSolverFactory ();

namespace Ifpack2 {
namespace Details {

void
registerLinearSolverFactory ()
{
  // Fill in the body of the function with all the type-specific
  // run-time registration functions.
  IFPACK2_INSTANTIATE_SLGN( LCL_CALL )
}

} // namespace Details
} // namespace Ifpack2


