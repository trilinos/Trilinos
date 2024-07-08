// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef AMESOS2_NEWMULTIVEC_ADAPTER_CPP
#define AMESOS2_NEWMULTIVEC_ADAPTER_CPP

#include "Amesos2_NewMultiVecAdapter_decl.hpp"

#ifdef HAVE_AMESOS2_EXPLICIT_INSTANTIATION
#  include "Amesos2_NewMultiVecAdapter_def.hpp"
#  include "Teuchos_ExplicitInstantiationHelpers.hpp"
namespace Amesos {
/* Need to figure out the explicit instantiation system for Amesos2, since we
 * are instantiating not only on the Scalar types, but Matrices and Vectors
 */

// TEUCHOS_CLASS_TEMPLATE_INSTANT_SCALAR_TYPES( )
}
#endif  // HAVE_AMESOS2_EXPLICIT_INSTANTIATION

#endif  // AMESOS2_NEWMULTIVEC_ADAPTER_CPP
