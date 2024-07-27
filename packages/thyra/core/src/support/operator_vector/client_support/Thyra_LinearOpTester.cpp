// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_LinearOpTester_decl.hpp"

#ifdef HAVE_THYRA_EXPLICIT_INSTANTIATION

#include "Thyra_LinearOpTester_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace Thyra {

TEUCHOS_CLASS_TEMPLATE_INSTANT_SCALAR_TYPES(LinearOpTester)

} // namespace Thyra

#endif // HAVE_THYRA_EXPLICIT_INSTANTIATION
