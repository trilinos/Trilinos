// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "TpetraCore_config.h"

#if defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)

// We protect the contents of this file with macros, to assist
// applications that circumvent Trilinos' build system.  (We do NOT
// recommend this.)  That way, they can still build this file, but as
// long as the macros have correct definitions, they won't build
// anything that's not enabled.

#include "Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp"
#include "Tpetra_Details_getGraphDiagOffsets_decl.hpp"
#include "Tpetra_Details_getGraphDiagOffsets_def.hpp"
#include "TpetraCore_ETIHelperMacros.h"

namespace Tpetra {

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN( TPETRA_DETAILS_IMPL_GETGRAPHDIAGOFFSETS_INSTANT )

} // namespace Tpetra

#endif // Whether we should build this specialization
