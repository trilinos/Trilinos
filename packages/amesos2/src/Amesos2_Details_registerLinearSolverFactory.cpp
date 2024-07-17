// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Amesos2_Details_registerLinearSolverFactory.hpp"
#include "Amesos2_Details_LinearSolverFactory.hpp"
// Amesos2 has a required dependency on Tpetra, so we don't need to
// protect inclusion of Tpetra header files with a macro.
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Operator.hpp"

#ifdef HAVE_AMESOS2_EPETRA
#  include "Epetra_MultiVector.h"
#  include "Epetra_Operator.h"
#endif // HAVE_AMESOS2_EPETRA

#include "TpetraCore_ETIHelperMacros.h"

// Define Tpetra instantiation macros and typedefs that make the
// macros work.  The fix for Bug 6380 makes this work whether or not
// ETI is ON.  We use the Tpetra macros because Amesos2 doesn't have
// its own macros.
TPETRA_ETI_MANGLING_TYPEDEFS()

// Macro that registers Amesos2's LinearSolverFactory for Tpetra
// objects, for the given four template parameters (Scalar = SC,
// LocalOrdinal = LO, GlobalOrdinal = GO, Node = NT).  The macro is
// local to this file.
//
// NOTE: This macro does NOT do explicit instantiation!  That's why I
// call it LCL_CALL and not LCL_INST.  We are just using the macros to
// invoke this class method over the set of enabled template
// parameters.
#define LCL_CALL( SC, LO, GO, NT ) \
  ::Amesos2::Details::LinearSolverFactory< ::Tpetra::MultiVector<SC, LO, GO, NT>, \
                                           ::Tpetra::Operator<SC, LO, GO, NT>, \
                                           typename ::Tpetra::MultiVector<SC, LO, GO, NT>::mag_type>::registerLinearSolverFactory ();

namespace Amesos2 {
namespace Details {

void
registerLinearSolverFactory ()
{
  // Fill in the body of the function with all the type-specific
  // run-time registration functions, for registering Amesos2's
  // LinearSolverFactory with Tpetra objects.
  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( LCL_CALL )

  // If Epetra is enabled in Amesos2, also register Amesos2's
  // LinearSolverFactory for Epetra objects.
#ifdef HAVE_AMESOS2_EPETRA
  ::Amesos2::Details::LinearSolverFactory<Epetra_MultiVector,
    Epetra_Operator, double>::registerLinearSolverFactory ();
#endif // HAVE_AMESOS2_EPETRA
}

} // namespace Details
} // namespace Amesos2


