// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef AMESOS2_TYPEMAP_HPP
#define AMESOS2_TYPEMAP_HPP

#include <Teuchos_ScalarTraits.hpp>

namespace Amesos2 {

/** \brief Map types to solver-specific data-types and enums.
 * 
 * \struct TypeMap
 *
 * Some direct linear sparse solvers have custom data types that are
 * more commonly represented as other data types.  For example,
 * Superlu uses a custom \c doublecomplex data types to represent
 * double-precision complex data.  Such a scalar is more commonly
 * represented by \c std::complex<double> .  Superlu then also uses an
 * \c enum class called \c Dtype to flag the data type for certain
 * methods.  Amesos2 uses TypeMap to easily access such data types.
 *
 * This class can be template specialized for each Amesos2::Solver
 * subclass and its supported types.  This also provides a compile
 * time check for whether a solver can handle a data type.  It is up
 * to the programmer/user to determine whether the data-type they have
 * may be safely coerced to a type supported by the ConcreteSolver,
 * perhaps with the help of Teuchos::as<>.  Appropriate type
 * conversions may be provided with through template specialization of
 * Teuchos::as<>.
 *
 * The default instance is empty, but specialized instances for each
 * ConcreteSolver should contain at the minimum a \c typedef called \c
 * type and other typedefs as appropriate for the ConcreteSolver's
 * needs
 * 
 * \tparam ConcreteSolver A Amesos2::Solver type for which these mappings hold
 * \tparam Scalar The Scalar type that is being mapped
 */
template <template<class,class> class ConcreteSolver,
          typename Scalar>
struct TypeMap {};


} // end namespace Amesos2

#endif  // AMESOS2_TYPEMAP_HPP
