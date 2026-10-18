// @HEADER
// *****************************************************************************
//      TpetraExt: Tpetra Extended - Linear Algebra Services Package
//
// Copyright 2025 NTESS
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "TpetraExt_PointToBlockDiagPermute_decl.hpp"

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION

#include "TpetraCore_ETIHelperMacros.h"
#include "TpetraExt_PointToBlockDiagPermute_def.hpp"

namespace Tpetra::Ext {

TPETRA_ETI_MANGLING_TYPEDEFS()

#define TPETRA_POINTTOBLOCKDIAGPERMUTE_INSTANT(SCALAR, LO, GO, NODE) \
  template class PointToBlockDiagPermute<SCALAR, LO, GO, NODE>;

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(TPETRA_POINTTOBLOCKDIAGPERMUTE_INSTANT)

}  // namespace Tpetra::Ext

#endif  // HAVE_TPETRA_EXPLICIT_INSTANTIATION
