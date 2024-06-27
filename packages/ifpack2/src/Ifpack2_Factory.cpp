// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Factory_decl.hpp"

#ifdef HAVE_IFPACK2_EXPLICIT_INSTANTIATION
#  include "Ifpack2_Factory_def.hpp"
#  include "Ifpack2_ExplicitInstantiationHelpers.hpp"
#  include "Ifpack2_ETIHelperMacros.h"
#endif // HAVE_IFPACK2_EXPLICIT_INSTANTIATION

namespace Ifpack2 {


bool supportsUnsymmetric (const std::string& prec_type)
{
  bool result = false;
  if (prec_type == "RELAXATION" ||
      prec_type == "CHEBYSHEV"  ||
      prec_type == "DIAGONAL"   ||
      prec_type == "MDF"        ||
      prec_type == "RILUK"      ||
      prec_type == "RBILUK"     ||
      prec_type == "ILUT"       ||
      prec_type == "SCHWARZ"    ||
      prec_type == "KRYLOV"
#ifdef HAVE_IFPACK_HYPRE
      || prec_type == "HYPRE"
#endif
      )
  {
    result = true;
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::invalid_argument, "Ifpack2::supportsUnsymmetric: "
      "Unrecognized preconditioner type prec_type = \"" << prec_type
      << "\"");
  }
  return result;
}

#ifdef HAVE_IFPACK2_EXPLICIT_INSTANTIATION

// We can't use the usual IFPACK2_* class macro here because
// OneLevelFactory is not a templated class; its methods are.
#define LCLINST(S, LO, GO)

  IFPACK2_ETI_MANGLING_TYPEDEFS()

  IFPACK2_INSTANTIATE_SLG_REAL( LCLINST )

#endif // HAVE_IFPACK2_EXPLICIT_INSTANTIATION

} // namespace Ifpack2

