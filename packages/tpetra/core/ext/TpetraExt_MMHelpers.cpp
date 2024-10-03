/*
// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

*/

#include "TpetraExt_MMHelpers.hpp"

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION

#include "TpetraCore_ETIHelperMacros.h"
#include "TpetraExt_MMHelpers_def.hpp"

namespace Tpetra {

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN(TPETRA_CRSMATRIXSTRUCT_INSTANT)

  TPETRA_INSTANTIATE_SLGN(TPETRA_CRSWRAPPER_INSTANT)

  TPETRA_INSTANTIATE_SLGN(TPETRA_CRSWRAPPER_CRSMATRIX_INSTANT)

  TPETRA_INSTANTIATE_SLGN(TPETRA_CRSWRAPPER_GRAPHBUILDER_INSTANT)

} // namespace Tpetra

#endif // HAVE_TPETRA_EXPLICIT_INSTANTIATION
