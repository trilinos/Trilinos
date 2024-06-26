// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "BelosConfigDefs.hpp"

#ifdef HAVE_BELOS_TPETRA
#  include "BelosTpetraAdapter.hpp" // must go first
#  include "Belos_Details_registerLinearSolverFactory.hpp"
#  include "Belos_Details_LinearSolverFactory.hpp"
#  include "TpetraCore_ETIHelperMacros.h"

// Define Tpetra instantiation macros and typedefs that make the
// macros work.  The fix for Bug 6380 makes this work whether or not
// ETI is ON.  We use the Tpetra macros because Belos doesn't have
// its own macros.
TPETRA_ETI_MANGLING_TYPEDEFS()

// Macro that registers Belos's LinearSolverFactory for Tpetra
// objects, for the given four template parameters (Scalar = SC,
// LocalOrdinal = LO, GlobalOrdinal = GO, Node = NT).  The macro is
// local to this file.
//
// NOTE: This macro does NOT do explicit instantiation!  That's why I
// call it LCL_CALL and not LCL_INST.  We are just using the macros to
// invoke this class method over the set of enabled template
// parameters.
#define LCL_CALL( SC, LO, GO, NT ) \
  ::Belos::Details::LinearSolverFactory< ::Tpetra::MultiVector<SC, LO, GO, NT>, \
                                         ::Tpetra::Operator<SC, LO, GO, NT>, \
                                         SC, \
                                         typename ::Tpetra::MultiVector<SC, LO, GO, NT>::mag_type>::registerLinearSolverFactory ();

namespace Belos {
namespace Details {
namespace Tpetra {

void
registerLinearSolverFactory ()
{
  // Fill in the body of the function with all the type-specific
  // run-time registration functions, for registering Belos's
  // LinearSolverFactory with Tpetra objects.
  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( LCL_CALL )
}

} // namespace Tpetra
} // namespace Details
} // namespace Belos

#endif // HAVE_BELOS_TPETRA

