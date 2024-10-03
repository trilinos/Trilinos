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

#include "Tpetra_Map.hpp"

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION

#include "TpetraCore_ETIHelperMacros.h"
#include "Tpetra_Map_def.hpp"

namespace Tpetra {

  TPETRA_ETI_MANGLING_TYPEDEFS()

  // for all nodes, lo, go
  TPETRA_INSTANTIATE_LGN(TPETRA_MAP_INSTANT)

  // for default node, all lo,go
  TPETRA_INSTANTIATE_LG(TPETRA_MAP_INSTANT_DEFAULTNODE)

} // namespace Tpetra

#endif // HAVE_TPETRA_EXPLICIT_INSTANTIATION
