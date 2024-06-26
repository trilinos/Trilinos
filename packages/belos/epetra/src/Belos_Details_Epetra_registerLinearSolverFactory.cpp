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

#ifdef HAVE_BELOS_EPETRA
#  include "Epetra_MultiVector.h"
#  include "Epetra_Operator.h"
#  include "BelosEpetraAdapter.hpp"

namespace Belos {
namespace Details {
namespace Epetra {

void
registerLinearSolverFactory ()
{
  // If Epetra is enabled in Belos, also register Belos's
  // LinearSolverFactory for Epetra objects.
  ::Belos::Details::LinearSolverFactory<Epetra_MultiVector,
    Epetra_Operator, double, double>::registerLinearSolverFactory ();
}

} // namespace Epetra
} // namespace Details
} // namespace Belos

#endif // HAVE_BELOS_EPETRA
