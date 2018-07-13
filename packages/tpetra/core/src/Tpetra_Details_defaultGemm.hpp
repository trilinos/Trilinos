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

#ifndef TPETRA_DETAILS_DEFAULTGEMM_HPP
#define TPETRA_DETAILS_DEFAULTGEMM_HPP

/// \file Tpetra_Details_defaultGemm.hpp
/// \brief Default implementation of local (but process-global) GEMM
///   (dense matrix-matrix multiply), for Tpetra::MultiVector.
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

namespace { // (anonymous)

  template<class ValueType>
  ValueType
  conditionalConjugate (const ValueType& x, const bool conj)
  {
    return conj ? Kokkos::Details::ArithTraits<ValueType>::conj (x) : x;
  }

} // namespace (anonymous)

/// \brief Default implementation of dense matrix-matrix multiply on a
///   single MPI process: <tt>C := alpha*A*B + beta*C</tt>.
///
/// \tparam ViewType1 Type of the first matrix input A.
/// \tparam ViewType2 Type of the second matrix input B.
/// \tparam ViewType3 Type of the third matrix input/output C.
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
gemm (const char transA,
      const char transB,
      const CoefficientType& alpha,
      const ViewType1& A,
      const ViewType2& B,
      const CoefficientType& beta,
      const ViewType3& C)
{
  // Assert that A, B, and C are in fact matrices
  static_assert (ViewType1::rank == 2, "GEMM: A must have rank 2 (be a matrix).");
  static_assert (ViewType2::rank == 2, "GEMM: B must have rank 2 (be a matrix).");
  static_assert (ViewType3::rank == 2, "GEMM: C must have rank 2 (be a matrix).");

  typedef typename ViewType3::non_const_value_type c_value_type;
  typedef Kokkos::Details::ArithTraits<CoefficientType> STS;
  const CoefficientType ZERO = STS::zero ();
  const CoefficientType ONE = STS::one ();

  // Get the dimensions
  const IndexType m = C.extent (0);
  const IndexType n = C.extent (1);
  const IndexType k = (transA == 'N' || transA == 'n') ?
    A.extent (1) : A.extent (0);

  const bool conjA = transA == 'C' || transA == 'c';
  const bool conjB = transB == 'C' || transB == 'c';

  // quick return if possible
  if (alpha == ZERO && beta == ONE) {
    return;
  }

  // And if alpha equals zero...
  if (alpha == ZERO) {
    if (beta == ZERO) {
      for (IndexType i = 0; i < m; ++i) {
        for (IndexType j = 0; j < n; ++j) {
          C(i,j) = ZERO;
        }
      }
    }
    else {
      for (IndexType i = 0; i < m; ++i) {
        for (IndexType j = 0; j < n; ++j) {
          C(i,j) = beta*C(i,j);
        }
      }
    }
    return;
  }

  // Start the operations
  if (transB == 'n' || transB == 'N') {
    if (transA == 'n' || transA == 'N') {
      // Form C = alpha*A*B + beta*C
      for (IndexType j = 0; j < n; ++j) {
        if (beta == ZERO) {
          for (IndexType i = 0; i < m; ++i) {
            C(i,j) = ZERO;
          }
        }
        else if (beta != ONE) {
          for (IndexType i = 0; i < m; ++i) {
            C(i,j) = beta*C(i,j);
          }
        }
        for (IndexType l = 0; l < k; ++l) {
          // Don't use c_value_type here, since it unnecessarily
          // forces type conversion before we assign to C(i,j).
          auto temp = alpha * conditionalConjugate (B(l,j), conjB);
          for (IndexType i = 0; i < m; ++i) {
            C(i,j) = C(i,j) + temp * conditionalConjugate (A(i,l), conjA);
          }
        }
      }
    }
    else {
      // Form C = alpha*A**T*B + beta*C
      for (IndexType j = 0; j < n; ++j) {
        for (IndexType i = 0; i < m; ++i) {
          c_value_type temp = ZERO;
          for (IndexType l = 0; l < k; ++l) {
            temp = temp + conditionalConjugate (A(l,i), conjA) *
              conditionalConjugate (B(l,j), conjB);
          }
          if (beta == ZERO) {
            C(i,j) = alpha*temp;
          }
          else {
            C(i,j) = alpha*temp + beta * C(i,j);
          }
        }
      }
    }
  }
  else {
    if (transA == 'n' || transA == 'N') {
      // Form C = alpha*A*B**T + beta*C
      for (IndexType j = 0; j < n; ++j) {
        if (beta == ZERO) {
          for (IndexType i = 0; i < m; ++i) {
            C(i,j) = ZERO;
          }
        }
        else if (beta != ONE) {
          for (IndexType i = 0; i < m; ++i) {
            C(i,j) = beta*C(i,j);
          }
        }
        for (IndexType l = 0; l < k; ++l) {
          // Don't use c_value_type here, since it unnecessarily
          // forces type conversion before we assign to C(i,j).
          auto temp = alpha * conditionalConjugate (B(j,l), conjB);
          for (IndexType i = 0; i < m; ++i) {
            C(i,j) = C(i,j) + temp * conditionalConjugate (A(i,l), conjA);
          }
        }
      }
    }
    else {
      // Form C = alpha*A**T*B**T + beta*C
      for (IndexType j = 0; j < n; ++j) {
        for (IndexType i = 0; i < m; ++i) {
          c_value_type temp = ZERO;
          for (IndexType l = 0; l < k; ++l) {
            temp = temp + conditionalConjugate (A(l,i), conjA) *
              conditionalConjugate (B(j,l), conjB);
          }
          if (beta == ZERO) {
            C(i,j) = alpha*temp;
          }
          else {
            C(i,j) = alpha*temp + beta*C(i,j);
          }
        }
      }
    }
  }
}

} // namespace Default
} // namespace Blas
} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_DEFAULTGEMM_HPP
