// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

///
/// \file Epetra_TsqrMessenger.hpp
/// \brief Function for wrapping Epetra_Comm in a communicator wrapper
///   that TSQR can use.
///
/// \note (mfh 27 Oct 2010) This file is in Tpetra (rather than
///   Epetra, where it would seem to belong) as a temporary fix.
///   Otherwise, Epetra would need an optional package dependency on
///   Teuchos and Kokkos, which would break third-party code linking
///   to the Epetra library.  Third-party code should use FIND_PACKAGE
///   on Trilinos to get the correct list of libraries against which
///   to link, but we make this easy temporary fix now so they have
///   time to fix their build systems later.
///

#ifndef EPETRA_TSQRMESSENGER_HPP
#define EPETRA_TSQRMESSENGER_HPP

#include <Tpetra_ConfigDefs.hpp>

#if !defined(TPETRA_ENABLE_DEPRECATED_CODE)
#error This file is deprecated due to Epetra removal and will be removed
#endif

#endif // EPETRA_TSQRMESSENGER_HPP

