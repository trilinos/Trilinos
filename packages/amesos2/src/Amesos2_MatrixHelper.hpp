// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef AMESOS2_MATRIXHELPER_HPP
#define AMESOS2_MATRIXHELPER_HPP

// #include "Amesos2_MatrixAdapter.hpp"
// #include "Amesos2_MultiVecAdapter.hpp"

namespace Amesos2 {

/**
 * \brief convert Matrices and MultiVectors into the appropriate
 * format for a third-party solver.
 *
 * The particular functions that must be implemented, and the
 * signatures in each will vary depending on the third-party solver's
 * needs.  This templated \c struct will just provide a central
 * location for functions which deal with Matrix and MultiVector
 * conversions.
 *
 * \tparam ConcreteSolver The Amesos2::Solver instance for which the
 * functions hold.
 */
template <template <typename,typename> class ConcreteSolver>
struct MatrixHelper 
{};

} // end namespace Amesos2

#endif  // AMESOS2_MATRIXHELPER_HPP
