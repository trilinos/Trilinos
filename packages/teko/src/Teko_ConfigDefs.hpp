// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_ConfigDefs_hpp__
#define __Teko_ConfigDefs_hpp__

#include "Teko_Config.h"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map_decl.hpp"

namespace Teko {

typedef double ST;
typedef Tpetra::Map<>::global_ordinal_type GO;
typedef Tpetra::Map<>::local_ordinal_type LO;
typedef Tpetra::Map<>::node_type NT;

}  // namespace Teko

#endif
