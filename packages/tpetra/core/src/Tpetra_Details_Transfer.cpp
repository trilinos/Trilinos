// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_Transfer_decl.hpp"

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION

#include "TpetraCore_ETIHelperMacros.h"
#include "Tpetra_Details_Transfer_def.hpp"

namespace Tpetra {
namespace Details {

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_LGN(TPETRA_DETAILS_TRANSFER_INSTANT)

} // namespace Details
} // namespace Tpetra

#endif // HAVE_TPETRA_EXPLICIT_INSTANTIATION
