// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include "Teuchos_ConfigDefs.hpp"
#include "RTOpPack_RTOpT_decl.hpp"

#if defined(HAVE_RTOP_EXPLICIT_INSTANTIATION) && defined(HAVE_TEUCHOS_LONG_DOUBLE)

#include "RTOpPack_RTOpT_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace RTOpPack {

TEUCHOS_CLASS_TEMPLATE_INSTANT_LONG_DOUBLE(RTOpT)

} // namespace RTOpPack

#endif // HAVE_TEUCHOS_EXCPLICIT_INSTANTIATION && HAVE_TEUCHOS_LONG_DOUBLE
