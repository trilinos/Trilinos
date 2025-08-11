// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_ETI_4ARG_MV_HPP
#define XPETRA_ETI_4ARG_MV_HPP

// The macro "XPETRA_ETI_GROUP" must be defined prior to including this file.

// We need to define these typedefs as it is not possible to properly expand
// macros with colons in them
#if defined(HAVE_XPETRA_TPETRA)
#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>
TPETRA_ETI_MANGLING_TYPEDEFS()
#endif

// Epetra = on, Tpetra = off

// Epetra = on, Tpetra = on

// Epetra = off, Tpetra = on
#if !defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_TPETRA)
TPETRA_INSTANTIATE_MULTIVECTOR(XPETRA_ETI_GROUP)
#endif

#endif  // ifndef XPETRA_ETI_4ARG_MV_HPP
