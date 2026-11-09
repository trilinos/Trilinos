// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef __IFPACK2_TEST_ADDITIVESCHWARZ_AMESOS2_SOLVE_HPP
#define __IFPACK2_TEST_ADDITIVESCHWARZ_AMESOS2_SOLVE_HPP

#include "TypedefsAndIncludes.hpp"

namespace Ifpack2 {
namespace Test {

bool solveTwice(const sparse_mat_type& A,
                const sparse_mat_type& A2,
                const multivector_type& b,
                const multivector_type& b2,
                const int overlapLevel,
                const int numBlocks,
                const int maxIters,
                const double tol,
                const bool reorder,
                const std::string& innerPrecondType,
                const bool expectOneIter);

}  // namespace Test
}  // namespace Ifpack2

#endif  // __IFPACK2_TEST_ADDITIVESCHWARZ_AMESOS2_SOLVE_HPP
