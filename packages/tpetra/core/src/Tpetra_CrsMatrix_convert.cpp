// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_ConfigDefs.hpp"

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION

#include "Tpetra_CrsMatrix_decl.hpp"
#include "TpetraCore_ETIHelperMacros.h"
#include "Tpetra_CrsMatrix_def.hpp"
#include "Tpetra_CrsGraph_def.hpp"

namespace Tpetra {

 TPETRA_ETI_MANGLING_TYPEDEFS()

 TPETRA_INSTANTIATE_CONVERT(TPETRA_CRSMATRIX_CONVERT_INSTANT)

} // namespace Tpetra

#endif // HAVE_TPETRA_EXPLICIT_INSTANTIATION
