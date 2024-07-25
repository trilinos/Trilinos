// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Belos_Details_registerLinearSolverFactory.hpp"
#include "Belos_Details_LinearSolverFactory.hpp"
#include "Teuchos_ConfigDefs.hpp" // for __attribute__((weak)) check

// FIXME (mfh 23 Aug 2015) Belos' main library is upstream from
// libraries where Belos' Epetra and Tpetra specializations live.
// That's why we need to use weak symbols here.  Otherwise, link
// errors result.  Fixing this requires the same run-time registration
// solution that Bug 6392 is all about.

#if defined(HAVE_BELOS_EPETRA) && defined(HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK)
namespace Belos {
namespace Details {
namespace Epetra {
  extern void __attribute__((weak)) registerLinearSolverFactory ();
} // namespace Epetra
} // namespace Details
} // namespace Belos
#endif // defined(HAVE_BELOS_EPETRA) && defined(HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK)

#if defined(HAVE_BELOS_TPETRA) && defined(HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK)
namespace Belos {
namespace Details {
namespace Tpetra {
  extern void __attribute__((weak)) registerLinearSolverFactory ();
} // namespace Tpetra
} // namespace Details
} // namespace Belos
#endif // defined(HAVE_BELOS_TPETRA) && defined(HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK)

//
// FIXME (mfh 23 Aug 2015) We should add Thyra as well.
//

namespace Belos {
namespace Details {

void
registerLinearSolverFactory ()
{
#if defined(HAVE_BELOS_TPETRA) && defined(HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK)
  // It's a weak symbol, so it might be NULL.
  if (::Belos::Details::Tpetra::registerLinearSolverFactory != NULL) {
    ::Belos::Details::Tpetra::registerLinearSolverFactory ();
  }
#endif // defined(HAVE_BELOS_TPETRA) && defined(HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK)

#if defined(HAVE_BELOS_EPETRA) && defined(HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK)
  // It's a weak symbol, so it might be NULL.
  if (::Belos::Details::Epetra::registerLinearSolverFactory != NULL) {
    ::Belos::Details::Epetra::registerLinearSolverFactory ();
  }
#endif // defined(HAVE_BELOS_EPETRA) && defined(HAVE_TEUCHOS_CXX_ATTRIBUTE_WEAK)
}

} // namespace Details
} // namespace Belos


