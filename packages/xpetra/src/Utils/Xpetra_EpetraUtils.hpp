// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_COMM_HPP
#define XPETRA_COMM_HPP

//! Conversion between Epetra and Teuchos objects

#include "Xpetra_ConfigDefs.hpp"

#ifdef HAVE_XPETRA_EPETRA

// header files for comm objects conversion
#include <Teuchos_Comm.hpp>
#include <Epetra_Comm.h>

// header file for Teuchos::ETransp
#include <Teuchos_BLAS_types.hpp>

namespace Xpetra {

using Teuchos::RCP;

//! Convert a Teuchos_Comm to an Epetra_Comm.
const RCP<const Epetra_Comm> toEpetra(const RCP<const Teuchos::Comm<int> >& comm);

//! Convert an Epetra_Comm.to a Teuchos_Comm.
const RCP<const Teuchos::Comm<int> > toXpetra(const Epetra_Comm& comm);

//! Convert a Teuchos::ETransp to an Epetra boolean.
bool toEpetra(Teuchos::ETransp);

}  // namespace Xpetra
#endif  // HAVE_XPETRA_EPETRA

#endif  // XPETRA_EPETRACOMM_HPP

// TODO: remove return RCP for toEpetra?
