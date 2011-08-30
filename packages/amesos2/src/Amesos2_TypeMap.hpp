// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package 
//                  Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
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
