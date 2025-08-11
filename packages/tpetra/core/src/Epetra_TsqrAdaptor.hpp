// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef EPETRA_TSQRADAPTOR_HPP
#define EPETRA_TSQRADAPTOR_HPP

/// \file Epetra_TsqrAdaptor.hpp
/// \brief Epetra_MultiVector to TSQR adaptor
///
/// \note (mfh 27 Oct 2010) This file is in Tpetra (rather than
/// Epetra, where it would seem to belong) as a temporary fix.
/// Otherwise, Epetra would need an optional package dependency on
/// Teuchos and Kokkos, which would break third-party code linking to
/// the Epetra library.  Third-party code should use FIND_PACKAGE on
/// Trilinos to get the correct list of libraries against which to
/// link, but we make this easy temporary fix now so they have time to
/// fix their build systems later.

#include "Tpetra_ConfigDefs.hpp"

#if defined(TPETRA_ENABLE_DEPRECATED_CODE)
#if defined(TPETRA_DEPRECATED_DECLARATIONS)
#warning This file is deprecated due to Epetra removal and will be removed
#endif
#else
#error This file is deprecated due to Epetra removal and will be removed
#endif

#endif // EPETRA_TSQRADAPTOR_HPP

