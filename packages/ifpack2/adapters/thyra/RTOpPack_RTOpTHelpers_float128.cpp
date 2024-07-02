// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include "Teuchos_ConfigDefs.hpp"
#include "RTOpPack_RTOpTHelpers_decl.hpp"

#if defined(HAVE_RTOP_EXPLICIT_INSTANTIATION) && defined(HAVE_TEUCHOSCORE_QUADMATH)

#include "RTOpPack_RTOpTHelpers_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace RTOpPack {

  TEUCHOS_MACRO_TEMPLATE_INSTANT_FLOAT128(RTOPPACK_RTOPT_HELPERS_INSTANT_SCALAR)

} // namespace RTOpPack

#endif // HAVE_RTOP_EXPLICIT_INSTANTIATION && HAVE_TEUCHOSCORE_QUADMATH
