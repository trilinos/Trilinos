// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_DefaultProductVector_decl.hpp"

#if (defined(HAVE_THYRA_EXPLICIT_INSTANTIATION) || defined(THYRA_DEFAULT_PRODUCT_VECTOR_EXPLICIT_INSTANTIATION))

#include "Thyra_DefaultProductVector_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace Thyra {

TEUCHOS_MACRO_TEMPLATE_INSTANT_SCALAR_TYPES(THYRA_DEFAULT_PRODUCT_VECTOR_INSTANT)

}  // namespace Thyra

#endif  // HAVE_THYRA_EXPLICIT_INSTANTIATION
