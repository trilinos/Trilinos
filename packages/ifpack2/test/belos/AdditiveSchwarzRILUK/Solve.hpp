// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef __IFPACK2_TEST_ADDITIVESCHWARZ_RILUK_SOLVE_HPP
#define __IFPACK2_TEST_ADDITIVESCHWARZ_RILUK_SOLVE_HPP

#include "TypedefsAndIncludes.hpp"

namespace Ifpack2 {
namespace Test {

void
solve (const sparse_mat_type& A,
       const multivector_type& b,
       multivector_type& x,
       const int ilukFillLevel,
       const int overlapLevel,
       const int numBlocks,
       const int maxIters,
       const double tol,
       const bool reorder,
       const std::string& innerPrecondType);

} // namespace Test
} // namespace Ifpack2

#endif // __IFPACK2_TEST_ADDITIVESCHWARZ_RILUK_SOLVE_HPP
