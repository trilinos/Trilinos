// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MueLu_Details_registerLinearSolverFactory.hpp"
#include "MueLu_Details_LinearSolverFactory.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Operator.hpp"
#ifdef HAVE_MUELU_EPETRA
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#endif
#include "TpetraCore_ETIHelperMacros.h"

// Define Tpetra instantiation macros and typedefs that make the
// macros work.  The fix for Bug 6380 makes this work whether or not
// ETI is ON.  We use the Tpetra macros because MueLu doesn't have
// its own macos.
TPETRA_ETI_MANGLING_TYPEDEFS()

// Macro that registers MueLu's LinearSolverFactory for Tpetra
// objects, for the given four template parameters (Scalar = SC,
// LocalOrdinal = LO, GlobalOrdinal = GO, Node = NT).  The macro is
// local to this file.
//
// NOTE: This macro does NOT do explicit instantiation!  That's why I
// call it LCL_CALL and not LCL_INST.  We are just using the macros to
// invoke this class method over the set of enabled template
// parameters.
#define LCL_CALL(SC, LO, GO, NT)                                             \
  ::MueLu::Details::LinearSolverFactory<Tpetra::MultiVector<SC, LO, GO, NT>, \
                                        Tpetra::Operator<SC, LO, GO, NT>,    \
                                        typename Tpetra::MultiVector<SC, LO, GO, NT>::mag_type>::registerLinearSolverFactory();

namespace MueLu {
namespace Details {

void registerLinearSolverFactory() {
  // Fill in the body of the function with all the type-specific
  // run-time registration functions, for registering MueLu's
  // LinearSolverFactory with Tpetra objects.
  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(LCL_CALL)

  // If Epetra is enabled in MueLu, also register MueLu's
  // LinearSolverFactory for Epetra objects.
#ifdef HAVE_MUELU_EPETRA
  ::MueLu::Details::LinearSolverFactory<Epetra_MultiVector,
                                        Epetra_Operator, double>::registerLinearSolverFactory();
#endif  // HAVE_MUELU_EPETRA
}

}  // namespace Details
}  // namespace MueLu