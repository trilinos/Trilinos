// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_DETAILS_TPETRA_REGISTERLINEARSOLVERFACTORY_HPP
#define BELOS_DETAILS_TPETRA_REGISTERLINEARSOLVERFACTORY_HPP

#include "BelosConfigDefs.hpp"

#ifdef HAVE_BELOS_TPETRA

namespace Belos {
namespace Details {
namespace Tpetra {

/// \brief Register Belos' LinearSolverFactory for Tpetra objects.
void registerLinearSolverFactory ();

} // namespace Tpetra
} // namespace Details
} // namespace Belos

#endif // HAVE_BELOS_TPETRA

#endif // BELOS_DETAILS_TPETRA_REGISTERLINEARSOLVERFACTORY_HPP
