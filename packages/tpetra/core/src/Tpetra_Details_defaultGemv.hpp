/*
//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// ************************************************************************
//@HEADER
*/

#ifndef TPETRA_DETAILS_DEFAULTGEMV_HPP
#define TPETRA_DETAILS_DEFAULTGEMV_HPP

/// \file Tpetra_Details_defaultGemv.hpp
/// \brief Default implementation of local (but process-global) GEMV
///   (dense matrix-vector multiply), for Tpetra::MultiVector.
///
/// \warning This file, and its contents, are an implementation detail
///   of Tpetra::MultiVector.  Either may disappear or change at any
///   time.

#include "TpetraCore_config.h"
#include "Kokkos_ArithTraits.hpp"
#include "Kokkos_Complex.hpp"
#include <type_traits>

namespace Tpetra {
namespace Details {
namespace Blas {
namespace Default {

/// \brief Default implementation of dense matrix-vector multiply on a
///   single MPI process: <tt>y := alpha*A*x + beta*y</tt>.
///
/// \tparam ViewType1 Type of the first matrix input A.
/// \tparam ViewType2 Type of the vector input x.
/// \tparam ViewType3 Type of the vector input/output y.
/// \tparam CoefficientType Type of the scalar coefficients alpha and beta.
/// \tparam IndexType Type of the index used in for loops; defaults to \c int.
///
/// ViewType1, ViewType2, and ViewType3 are all Kokkos::View specializations.
template<class ViewType1,
         class ViewType2,
         class ViewType3,
         class CoefficientType,
         class IndexType = int>
void
gemv (const char trans,
      const CoefficientType& alpha,
      const ViewType1& A,
      const ViewType2& x,
      const CoefficientType& beta,
      const ViewType3& y)
{
  // Assert that A, B, and C are in fact matrices
  static_assert (ViewType1::rank == 2, "GEMV: A must have rank 2 (be a matrix).");
  static_assert (ViewType2::rank == 1, "GEMV: x must have rank 1 (be a vector).");
  static_assert (ViewType3::rank == 1, "GEMV: y must have rank 1 (be a vector).");

  using Kokkos::ArithTraits;
  typedef typename ViewType1::non_const_value_type A_value_type;
  typedef typename ViewType3::non_const_value_type y_value_type;
  const CoefficientType ZERO = ArithTraits<CoefficientType>::zero ();
  const CoefficientType ONE = ArithTraits<CoefficientType>::one ();

  // // Get the dimensions
  // IndexType m, n;
  // if (trans == 'N' || trans == 'n') {
  //   m = static_cast<IndexType> (A.extent (0));
  //   n = static_cast<IndexType> (A.extent (1));
  // }
  // else {
  //   m = static_cast<IndexType> (A.extent (1));
  //   n = static_cast<IndexType> (A.extent (0));
  // }

  // quick return if possible
  if (alpha == ZERO && beta == ONE) {
    return; // y := 1*y
  }

  const bool noTrans = trans == 'n' || trans == 'N';
  if (noTrans) {
    if (alpha == ZERO) {
      if (beta == ZERO) {
        Kokkos::deep_copy (y, ArithTraits<y_value_type>::zero ());
      }
      else { // alpha == 1, beta != 0
        for (IndexType i = 0; i < static_cast<IndexType> (y.extent (0)); ++i) {
          y(i) = beta * y(i);
        }
      }
    }
    else { // alpha != 0
      for (IndexType i = 0; i < static_cast<IndexType> (y.extent (0)); ++i) {
        y_value_type y_i = (beta == ZERO) ?
          ArithTraits<y_value_type>::zero () :
          static_cast<y_value_type> (beta * y(i));
        for (IndexType j = 0; j < static_cast<IndexType> (x.extent (0)); ++j) {
          y_i += alpha * A(i,j) * x(j);
        }
        y(i) = y_i;
      }
    }
  }
  else { // ! noTrans
    const bool conj = ! (trans == 't' || trans == 'T');

    if (alpha == ZERO) {
      if (beta == ZERO) {
        Kokkos::deep_copy (y, ArithTraits<y_value_type>::zero ());
      }
      else { // alpha == 1, beta != 0
        for (IndexType j = 0; j < static_cast<IndexType> (y.extent (0)); ++j) {
          y(j) = beta * y(j);
        }
      }
    }
    else { // alpha != 0
      for (IndexType j = 0; j < static_cast<IndexType> (y.extent (0)); ++j) {
        y_value_type y_j = (beta == ZERO) ?
          ArithTraits<y_value_type>::zero () :
          static_cast<y_value_type> (beta * y(j));
        for (IndexType i = 0; i < static_cast<IndexType> (x.extent (0)); ++i) {
          const auto A_ij = conj ? Kokkos::ArithTraits<A_value_type>::conj (A(i,j)) : A(i,j);
          y_j += alpha * A_ij * x(i);
        }
        y(j) = y_j;
      }
    }
  }
}

} // namespace Default
} // namespace Blas
} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_DEFAULTGEMV_HPP
