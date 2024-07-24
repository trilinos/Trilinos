// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "TpetraCore_config.h"

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION

#include "Tpetra_BlockCrsMatrix_Helpers_decl.hpp"
#include "Tpetra_BlockCrsMatrix_Helpers_def.hpp"

#include "TpetraCore_ETIHelperMacros.h"

namespace Tpetra {
  
  TPETRA_ETI_MANGLING_TYPEDEFS()
  
  TPETRA_INSTANTIATE_LGN( TPETRA_CREATEMESHMAP_INSTANT )
  
} // namespace Tpetra

#endif // HAVE_TPETRA_INSTANTIATION
