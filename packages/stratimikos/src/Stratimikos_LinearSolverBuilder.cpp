// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stratimikos_LinearSolverBuilder_decl.hpp"

#ifdef HAVE_STRATIMIKOS_EXPLICIT_INSTANTIATION

#include "Stratimikos_LinearSolverBuilder_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace Stratimikos {

TEUCHOS_MACRO_TEMPLATE_INSTANT_SCALAR_TYPES(STRATIMIKOS_LINEARSOLVERBUILDER_INSTANT)

int LinearSolverBuilderHelpers::existingNameIndex(
  const Teuchos::ArrayView<std::string> namesArray, const std::string &name)
{
  typedef Teuchos::ArrayView<std::string>::const_iterator iter_t;
  const iter_t iter_begin = namesArray.begin(), iter_end = namesArray.end();
  const iter_t iter = std::find(iter_begin, iter_end, name);
  if (iter != iter_end) {
    return (iter - iter_begin);
  }
  return -1;
}

} // namespace Stratimikos

#endif // HAVE_STRATIMIKOS_EXPLICIT_INSTANTIATION
