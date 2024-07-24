// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Tpetra_Details_initializeKokkos.hpp
/// \brief Declaration of Tpetra::Details::initializeKokkos
///
/// \warning This is an implementation detail of Tpetra.  This header
///   file and/or its contents may change or disappear at any time.

#ifndef TPETRA_DETAILS_INITIALIZEKOKKOS_HPP
#define TPETRA_DETAILS_INITIALIZEKOKKOS_HPP

namespace Tpetra {
namespace Details {

/// \brief Initialize Kokkos, using command-line arguments (if any)
///   given to Teuchos::GlobalMPISession.
///
/// \warning This is an implementation detail of Tpetra.  This header
///   file and/or its contents may change or disappear at any time.
///
/// \warning Prefer Tpetra::initialize.  Give it the actual
///   command-line arguments.
void
initializeKokkos ();

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_INITIALIZEKOKKOS_HPP
