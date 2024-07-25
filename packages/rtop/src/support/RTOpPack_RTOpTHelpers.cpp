// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#include "RTOpPack_RTOpTHelpers_decl.hpp"


#ifdef RTOPPACK_RTOPT_HELPER_DUMP_OUTPUT
bool RTOpPack::rtop_helpers_dump_all = false;
#endif


#ifdef HAVE_RTOP_EXPLICIT_INSTANTIATION


#include "RTOpPack_RTOpTHelpers_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"


namespace RTOpPack {


TEUCHOS_MACRO_TEMPLATE_INSTANT_SCALAR_TYPES(
  RTOPPACK_RTOPT_HELPERS_INSTANT_SCALAR)


RTOPPACK_RTOPT_HELPERS_DEFAULTREDUCTTARGET_INSTANT(Ordinal)


} // namespace RTOpPack


#endif // HAVE_TEUCHOS_EXCPLICIT_INSTANTIATION
