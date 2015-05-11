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

/// \file Kokkos_Sparse.hpp
/// \brief Public interface to local computational kernels on sparse
///   matrices.
///
/// KokkosSparse::spmv implements local sparse matrix-vector multiply.
/// It computes y = beta*y + alpha*Op(A)*x, where x and y are either
/// both rank 1 (single vectors) or rank 2 (multivectors) Kokkos::View
/// instances, A is a KokkosSparse::CrsMatrix, and Op(A) is determined
/// by the \c mode input (either no transpose, transpose, or conjugate
/// transpose).  If beta == 0, ignore and overwrite the initial
/// entries of y; if alpha == 0, ignore the entries of A and x.
///
/// KokkosSparse::trsv implements local sparse triangular solve.
/// It solves Ax=b, where A is either upper or lower triangular.

#ifndef KOKKOS_SPARSE_HPP_
#define KOKKOS_SPARSE_HPP_

#ifdef KOKKOS_HAVE_CXX11
#include <type_traits>
#endif // KOKKOS_HAVE_CXX11

#include <Kokkos_Sparse_CrsMatrix.hpp>
#include <Kokkos_Sparse_impl_spmv.hpp>
#include <Kokkos_Sparse_trsv.hpp>

namespace KokkosSparse {

namespace {
  struct RANK_ONE{};
  struct RANK_TWO{};
}

template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
void
spmv(const char mode[],
     const AlphaType& alpha,
     const AMatrix& A,
     const XVector& x,
     const BetaType& beta,
     const YVector& y,
     const RANK_ONE)
{

#ifdef KOKKOS_HAVE_CXX11
  // Make sure that both x and y have the same rank.
  static_assert ((int) XVector::rank == (int) YVector::rank,
                 "KokkosBlas::spmv: Vector ranks do not match.");
  // Make sure that x (and therefore y) is rank 1.
  static_assert ((int) XVector::rank == 1, "KokkosBlas::spmv: "
                 "Both Vector inputs must have rank 1 or 2.");
  // Make sure that y is non-const.
  static_assert (Kokkos::Impl::is_same<typename YVector::value_type,
                                       typename YVector::non_const_value_type>::value,
                 "KokkosBlas::spmv: Output Vector must be non-const.");

#else
  // We prefer to use C++11 static_assert, because it doesn't give
  // "unused typedef" warnings, like the constructs below do.
  //
  // Make sure that both x and y have the same rank.
  typedef typename
    Kokkos::Impl::StaticAssert<(int) XVector::rank == (int) YVector::rank>::type
    Blas1_spmv_vector_ranks_do_not_match;
  // Make sure that x (and therefore y) is rank 1.
  typedef typename
    Kokkos::Impl::StaticAssert<(int) XVector::rank == 1 >::type
    Blas1_spmv_vector_rank_not_one;
#endif // KOKKOS_HAVE_CXX11

  // Check compatibility of dimensions at run time.
  if((mode[0]==NoTranspose[0])||(mode[0]==Conjugate[0])) {
    if ((x.dimension_1 () != y.dimension_1 ()) ||
        (static_cast<size_t> (A.numCols ()) > static_cast<size_t> (x.dimension_0 ())) ||
        (static_cast<size_t> (A.numRows ()) > static_cast<size_t> (y.dimension_0 ()))) {
      std::ostringstream os;
      os << "KokkosBlas::spmv: Dimensions do not match: "
         << ", A: " << A.numRows () << " x " << A.numCols()
         << ", x: " << x.dimension_0 () << " x " << x.dimension_1()
         << ", y: " << y.dimension_0 () << " x " << y.dimension_1()
         ;

      Kokkos::Impl::throw_runtime_exception (os.str ());
    }
  } else {
    if ((x.dimension_1 () != y.dimension_1 ()) ||
        (static_cast<size_t> (A.numCols ()) > static_cast<size_t> (y.dimension_0 ())) ||
        (static_cast<size_t> (A.numRows ()) > static_cast<size_t> (x.dimension_0()))) {
      std::ostringstream os;
      os << "KokkosBlas::spmv: Dimensions do not match (transpose): "
         << ", A: " << A.numRows () << " x " << A.numCols()
         << ", x: " << x.dimension_0 () << " x " << x.dimension_1()
         << ", y: " << y.dimension_0 () << " x " << y.dimension_1()
         ;

      Kokkos::Impl::throw_runtime_exception (os.str ());
    }
  }

  typedef KokkosSparse::CrsMatrix<typename AMatrix::const_value_type,
    typename AMatrix::non_const_ordinal_type,
    typename AMatrix::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename AMatrix::const_size_type> AMatrix_Internal;
  typedef Kokkos::View<typename XVector::const_value_type*,
    typename XVector::array_layout,
    typename XVector::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>,
    typename XVector::specialize> XVector_Internal;
  typedef Kokkos::View<typename YVector::non_const_value_type*,
    typename YVector::array_layout,
    typename YVector::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename YVector::specialize> YVector_Internal;
  AMatrix_Internal A_i = A;
  XVector_Internal x_i = x;
  YVector_Internal y_i = y;

  return Impl::SPMV<typename AMatrix_Internal::value_type,
             typename AMatrix_Internal::ordinal_type,
             typename AMatrix_Internal::device_type,
             typename AMatrix_Internal::memory_traits,
             typename AMatrix_Internal::size_type,
             typename XVector_Internal::value_type*,
             typename XVector_Internal::array_layout,
             typename XVector_Internal::device_type,
             typename XVector_Internal::memory_traits,
             typename XVector_Internal::specialize,
             typename YVector_Internal::value_type*,
             typename YVector_Internal::array_layout,
             typename YVector_Internal::device_type,
             typename YVector_Internal::memory_traits,
             typename YVector_Internal::specialize>::spmv(mode,alpha,A,x,beta,y);
}

}

#endif

#include<Kokkos_Sparse_MV.hpp>
