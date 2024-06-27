// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_DETAILS_EPETRA_REGISTERLINEARSOLVERFACTORY_HPP
#define BELOS_DETAILS_EPETRA_REGISTERLINEARSOLVERFACTORY_HPP

#include "BelosConfigDefs.hpp"

#ifdef HAVE_BELOS_EPETRA

namespace Belos {
namespace Details {
namespace Epetra {

/// \brief Register Belos' LinearSolverFactory for Epetra objects.
void registerLinearSolverFactory ();

} // namespace Epetra
} // namespace Details
} // namespace Belos

#endif // HAVE_BELOS_EPETRA

#endif // BELOS_DETAILS_EPETRA_REGISTERLINEARSOLVERFACTORY_HPP
