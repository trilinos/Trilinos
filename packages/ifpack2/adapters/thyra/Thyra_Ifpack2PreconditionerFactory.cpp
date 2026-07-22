// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_Ifpack2PreconditionerFactory_decl.hpp"

#include "Ifpack2_ConfigDefs.hpp"

#ifdef HAVE_IFPACK2_EXPLICIT_INSTANTIATION

#include "Thyra_Ifpack2PreconditionerFactory_def.hpp"

#include "Ifpack2_ExplicitInstantiationHelpers.hpp"
#include "Ifpack2_ETIHelperMacros.h"

namespace Thyra {

#define LCLINST(S, LO, GO, NO) \
  IFPACK2_INST(Ifpack2PreconditionerFactory, S, LO, GO, NO)

IFPACK2_ETI_MANGLING_TYPEDEFS()

IFPACK2_INSTANTIATE_SLGN(LCLINST)

}  // namespace Thyra

#endif
