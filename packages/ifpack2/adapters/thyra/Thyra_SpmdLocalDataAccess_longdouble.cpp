// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include "Teuchos_ConfigDefs.hpp"
#include "Thyra_ConfigDefs.hpp"

#if defined(HAVE_THYRA_EXPLICIT_INSTANTIATION) && defined(HAVE_TEUCHOS_LONG_DOUBLE)

#include "Thyra_SpmdLocalDataAccess_decl.hpp"
#include "Thyra_SpmdLocalDataAccess_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace Thyra {

  TEUCHOS_MACRO_TEMPLATE_INSTANT_LONG_DOUBLE(THYRA_SPMD_LOCAL_DATA_ACCESS_INSTANT)

} // namespace Thyra

#endif // HAVE_THYRA_EXPLICIT_INSTANTIATION && HAVE_TEUCHOS_LONG_DOUBLE
