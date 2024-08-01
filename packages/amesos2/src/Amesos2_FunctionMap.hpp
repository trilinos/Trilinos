// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
   \file   Amesos2_FunctionMap.hpp
   \author Eric Bavier <etbavier@sandia.gov>
   \date   Thu Jul 15 22:04:41 2010

   \brief  Declaration of Function mapping class for Amesos2.
*/


#ifndef AMESOS2_FUNCTIONMAP_HPP
#define AMESOS2_FUNCTIONMAP_HPP

namespace Amesos2 {

/**
 * \brief Passes functions to TPL functions based on type.
 *
 * Helper class which passes on function calls to the appropriate
 * Solver function based on the type of its scalar template argument.
 *
 * Some Solvers have solver and matrix builder functions defined based
 * on data type.  One function for complex, one for double precision
 * complex, another for \c float , and yet another for \c double.  To
 * work elegantly with the Amesos2::Solver interface we want to be
 * able to perform a single function call which is appropriate for the
 * scalar type of the Matrix and MultiVectors that we are working
 * with.  The \c FunctionMap class provides that capability.
 *
 * The class template is specialized for each Solver and each data
 * type that it supports, and errors are thrown for other data types.
 */
template <template <typename,typename> class ConcreteSolver, typename Scalar>
struct FunctionMap
{};


} // end namespace Amesos2

#endif  // AMESOS2_FUNCTIONMAP_HPP
