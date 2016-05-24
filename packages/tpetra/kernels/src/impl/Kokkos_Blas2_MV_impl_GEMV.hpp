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
#ifndef KOKKOS_BLAS2_MV_IMPL_GEMV_HPP_
#define KOKKOS_BLAS2_MV_IMPL_GEMV_HPP_

#include "TpetraKernels_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"

namespace KokkosBlas {
namespace Impl {

// Functor for a single-level parallel_for version of nontranspose
// GEMV.  The functor parallelizes over rows of the input matrix A.
template<class AViewType,
         class XViewType,
         class YViewType,
         const int alphaPreset, // 0 or 1 are specific values; -1 means general
         const int betaPreset, // 0 or 1 are specific values; -1 means general
         class AlphaCoeffType = typename YViewType::non_const_value_type,
         class BetaCoeffType = typename YViewType::non_const_value_type,
         class IndexType = typename AViewType::size_type>
struct SingleLevelNontransposeGEMV {
  SingleLevelNontransposeGEMV (const AlphaCoeffType& alpha,
                               const AViewType& A,
                               const XViewType& x,
                               const BetaCoeffType& beta,
                               const YViewType& y) :
    alpha_ (alpha), A_ (A), x_ (x), beta_ (beta), y_ (y)
  {
    static_assert (Kokkos::Impl::is_view<AViewType>::value,
                   "AViewType must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XViewType>::value,
                   "XViewType must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YViewType>::value,
                   "YViewType must be a Kokkos::View.");
    static_assert (static_cast<int> (AViewType::rank) == 2,
                   "AViewType must have rank 2.");
    static_assert (static_cast<int> (XViewType::rank) == 1,
                   "XViewType must have rank 1.");
    static_assert (static_cast<int> (YViewType::rank) == 1,
                   "YViewType must have rank 1.");
    static_assert (std::is_integral<IndexType>::value,
                   "IndexType must be an integer.");
    static_assert (alphaPreset == 0 || alphaPreset == 1 || alphaPreset == -1,
                   "Invalid alphaPreset value; valid values are 0, 1, and -1.");
    static_assert (betaPreset == 0 || betaPreset == 1 || betaPreset == -1,
                   "Invalid betaPreset value; valid values are 0, 1, and -1.");
  }

  // i is the current row of the input matrix A, and the current row
  // of the output vector y.
  KOKKOS_INLINE_FUNCTION void
  operator () (const IndexType& i) const
  {
    typedef typename std::decay<decltype (y_[i]) >::type y_value_type;

    y_value_type y_i;
    if (betaPreset == 0) {
      y_i = Kokkos::Details::ArithTraits<y_value_type>::zero ();
    }
    else if (betaPreset == 1) {
      y_i = y_[i];
    }
    else { // beta_ != 0 and beta != 1
      y_i = beta_ * y_[i];
    }

    const IndexType numCols = A_.dimension_1 ();
    if (alphaPreset == 0) {
      ; // do nothing
    }
    else if (alphaPreset == 1) {
      for (IndexType j = 0; j < numCols; ++j) {
        y_i += A_(i,j) * x_(j);
      }
    }
    else { // alpha_ != 0 and alpha_ != 1
      for (IndexType j = 0; j < numCols; ++j) {
        y_i += alpha_ * A_(i,j) * x_(j);
      }
    }

    y_[i] = y_i;
  }

private:
  AlphaCoeffType alpha_;
  typename AViewType::const_type A_;
  typename XViewType::const_type x_;
  BetaCoeffType beta_;
  YViewType y_;
};

// Functor for a single-level parallel_reduce version of (conjugate)
// transpose GEMV.  The functor parallelizes over rows of the input
// matrix A and the input vector x.  The output vector y is the
// reduction result.
//
// WARNING: NOT RECOMMENDED FOR CUDA.  Reduction result may have
// arbitrary length.  This is bad on CUDA because the CUDA
// implementation of Kokkos::parallel_reduce may use shared memory for
// intermediate results.
template<class AViewType,
         class XViewType,
         class YViewType,
         const bool conj,
         const int alphaPreset, // 0 or 1 are specific values; -1 means general
         const int betaPreset, // 0 or 1 are specific values; -1 means general
         class AlphaCoeffType = typename YViewType::non_const_value_type,
         class BetaCoeffType = typename YViewType::non_const_value_type,
         class IndexType = typename AViewType::size_type>
struct SingleLevelTransposeGEMV {
  typedef typename YViewType::non_const_value_type y_value_type;
  typedef y_value_type value_type[];
  IndexType value_count; // Kokkos needs this for reductions w/ array results

  SingleLevelTransposeGEMV (const AlphaCoeffType& alpha,
                            const AViewType& A,
                            const XViewType& x,
                            const BetaCoeffType& beta,
                            const YViewType& y) :
    value_count (A.dimension_1 ()), alpha_ (alpha),
    A_ (A), x_ (x), beta_ (beta), y_ (y)
  {
    static_assert (Kokkos::Impl::is_view<AViewType>::value,
                   "AViewType must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XViewType>::value,
                   "XViewType must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YViewType>::value,
                   "YViewType must be a Kokkos::View.");
    static_assert (static_cast<int> (AViewType::rank) == 2,
                   "AViewType must have rank 2.");
    static_assert (static_cast<int> (XViewType::rank) == 1,
                   "XViewType must have rank 1.");
    static_assert (static_cast<int> (YViewType::rank) == 1,
                   "YViewType must have rank 1.");
    static_assert (std::is_integral<IndexType>::value,
                   "IndexType must be an integer.");
    static_assert (alphaPreset == 0 || alphaPreset == 1 || alphaPreset == -1,
                   "Invalid alphaPreset value; valid values are 0, 1, and -1.");
    static_assert (betaPreset == 0 || betaPreset == 1 || betaPreset == -1,
                   "Invalid betaPreset value; valid values are 0, 1, and -1.");
  }

public:
  KOKKOS_INLINE_FUNCTION void
  init (value_type y_cur) const
  {
    for (IndexType j = 0; j < value_count; ++j) {
      y_cur[j] = Kokkos::Details::ArithTraits<y_value_type>::zero ();
    }
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type dst,
        const volatile value_type src) const
  {
    for (IndexType j = 0; j < value_count; ++j) {
      dst[j] += src[j];
    }
  }

  KOKKOS_INLINE_FUNCTION void
  final (value_type y_result) const
  {
    using Kokkos::Details::ArithTraits;

    for (IndexType j = 0; j < value_count; ++j) {
      // Sum into initial y_ values; use beta as a pre-multiplier if nonzero.
      const y_value_type y_j =
        beta_ == ArithTraits<BetaCoeffType>::zero () ?
        ArithTraits<y_value_type>::zero () :
        beta_ * y_[j];
      y_[j] = y_j + y_result[j];
    }
  }

  KOKKOS_INLINE_FUNCTION void
  operator () (const IndexType& i, value_type y_cur) const
  {
    using Kokkos::Details::ArithTraits;
    typedef ArithTraits<typename AViewType::non_const_value_type> KAT;

    const auto x_i = x_(i);
    if (alphaPreset == 1) {
      for (IndexType j = 0; j < value_count; ++j) {
        const auto A_ij = conj ? KAT::conj (A_(i,j)) : A_(i,j);
        y_cur[j] += A_ij * x_i;
      }
    }
    else if (alphaPreset != 0) { // alpha_ != 0 and alpha_ != 1
      for (IndexType j = 0; j < value_count; ++j) {
        const auto A_ij = conj ? KAT::conj (A_(i,j)) : A_(i,j);
        y_cur[j] += alpha_ * A_ij * x_i;
      }
    }
  }

private:
  AlphaCoeffType alpha_;
  typename AViewType::const_type A_;
  typename XViewType::const_type x_;
  BetaCoeffType beta_;
  YViewType y_;
};

// Single-level parallel version of GEMV.
template<class AViewType,
         class XViewType,
         class YViewType,
         class AlphaCoeffType = typename YViewType::non_const_value_type,
         class BetaCoeffType = typename YViewType::non_const_value_type,
         class IndexType = typename AViewType::size_type>
void
singleLevelGemv (const char trans[],
                 const AlphaCoeffType& alpha,
                 const AViewType& A,
                 const XViewType& x,
                 const BetaCoeffType& beta,
                 const YViewType& y)
{
  static_assert (Kokkos::Impl::is_view<AViewType>::value,
                 "AViewType must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XViewType>::value,
                 "XViewType must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<YViewType>::value,
                 "YViewType must be a Kokkos::View.");
  static_assert (static_cast<int> (AViewType::rank) == 2,
                 "AViewType must have rank 2.");
  static_assert (static_cast<int> (XViewType::rank) == 1,
                 "XViewType must have rank 1.");
  static_assert (static_cast<int> (YViewType::rank) == 1,
                 "YViewType must have rank 1.");
  static_assert (std::is_integral<IndexType>::value,
                 "IndexType must be an integer");

  typedef typename YViewType::non_const_value_type y_value_type;
  typedef typename AViewType::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, IndexType> policy_type;

  policy_type range (0, A.dimension_0 ());
  const char tr = trans[0];

  // The transpose and conjugate transpose cases where A has zero rows
  // need special handling.  These are equivalent to y := beta*y.  We
  // could implement this using KokkosBlas::scal, but we don't want to
  // depend on that or its implementation details.  Instead, we reuse
  // an instantiation of the non-transpose case for alpha=0.
  if (A.dimension_0 () == 0 && (tr != 'N' && tr != 'n')) {
    if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::zero ()) {
      Kokkos::deep_copy (y, Kokkos::Details::ArithTraits<BetaCoeffType>::zero ());
    }
    else if (beta != Kokkos::Details::ArithTraits<BetaCoeffType>::one ()) {
      // "Fake out" a scal() by using the non-transpose alpha=0,
      // general beta case.  This assumes that the functor doesn't
      // check dimensions.
      typedef SingleLevelNontransposeGEMV<AViewType, XViewType, YViewType,
        0, -1, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
      functor_type functor (alpha, A, x, beta, y);
      Kokkos::parallel_for (policy_type (0, A.dimension_1 ()), functor);
    }
    return;
  }

  if (tr == 'N' || tr == 'n') {
    if (alpha == Kokkos::Details::ArithTraits<AlphaCoeffType>::zero ()) {
      if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::zero ()) {
        // Fill y with zeros
        Kokkos::deep_copy (y, Kokkos::Details::ArithTraits<y_value_type>::zero ());
      }
      else if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::one ()) {
        // Do nothing (y := 1 * y)
      }
      else { // beta != 0 && beta != 1
        typedef SingleLevelNontransposeGEMV<AViewType, XViewType, YViewType,
          0, -1, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
        functor_type functor (alpha, A, x, beta, y);
        Kokkos::parallel_for (range, functor);
      }
    }
    else if (alpha == Kokkos::Details::ArithTraits<AlphaCoeffType>::one ()) {
      if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::zero ()) {
        typedef SingleLevelNontransposeGEMV<AViewType, XViewType, YViewType,
          1, 0, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
        functor_type functor (alpha, A, x, beta, y);
        Kokkos::parallel_for (range, functor);
      }
      else if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::one ()) {
        typedef SingleLevelNontransposeGEMV<AViewType, XViewType, YViewType,
          1, 1, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
        functor_type functor (alpha, A, x, beta, y);
        Kokkos::parallel_for (range, functor);
      }
      else { // beta != 0 && beta != 1
        typedef SingleLevelNontransposeGEMV<AViewType, XViewType, YViewType,
          1, -1, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
        functor_type functor (alpha, A, x, beta, y);
        Kokkos::parallel_for (range, functor);
      }
    }
    else { // alpha != 0 and alpha != 1
      if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::zero ()) {
        typedef SingleLevelNontransposeGEMV<AViewType, XViewType, YViewType,
          -1, 0, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
        functor_type functor (alpha, A, x, beta, y);
        Kokkos::parallel_for (range, functor);
      }
      else if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::one ()) {
        typedef SingleLevelNontransposeGEMV<AViewType, XViewType, YViewType,
          -1, 1, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
        functor_type functor (alpha, A, x, beta, y);
        Kokkos::parallel_for (range, functor);
      }
      else { // beta != 0 && beta != 1
        typedef SingleLevelNontransposeGEMV<AViewType, XViewType, YViewType,
          -1, -1, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
        functor_type functor (alpha, A, x, beta, y);
        Kokkos::parallel_for (range, functor);
      }
    }
  }
  else if (tr == 'T' || tr == 't') { // transpose, no conjugate
    if (alpha == Kokkos::Details::ArithTraits<AlphaCoeffType>::zero ()) {
      if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::zero ()) {
        // Fill y with zeros
        Kokkos::deep_copy (y, Kokkos::Details::ArithTraits<y_value_type>::zero ());
      }
      else if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::one ()) {
        // Do nothing (y := 1 * y)
      }
      else { // beta != 0 && beta != 1
        typedef SingleLevelTransposeGEMV<AViewType, XViewType, YViewType,
          false, 0, -1, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
        functor_type functor (alpha, A, x, beta, y);
        Kokkos::parallel_reduce (range, functor);
      }
    }
    else if (alpha == Kokkos::Details::ArithTraits<AlphaCoeffType>::one ()) {
      if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::zero ()) {
        typedef SingleLevelTransposeGEMV<AViewType, XViewType, YViewType,
          false, 1, 0, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
        functor_type functor (alpha, A, x, beta, y);
        Kokkos::parallel_reduce (range, functor);
      }
      else if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::one ()) {
        typedef SingleLevelTransposeGEMV<AViewType, XViewType, YViewType,
          false, 1, 1, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
        functor_type functor (alpha, A, x, beta, y);
        Kokkos::parallel_reduce (range, functor);
      }
      else { // beta != 0 && beta != 1
        typedef SingleLevelTransposeGEMV<AViewType, XViewType, YViewType,
          false, 1, -1, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
        functor_type functor (alpha, A, x, beta, y);
        Kokkos::parallel_reduce (range, functor);
      }
    }
    else { // alpha != 0 and alpha != 1
      if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::zero ()) {
        typedef SingleLevelTransposeGEMV<AViewType, XViewType, YViewType,
          false, -1, 0, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
        functor_type functor (alpha, A, x, beta, y);
        Kokkos::parallel_reduce (range, functor);
      }
      else if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::one ()) {
        typedef SingleLevelTransposeGEMV<AViewType, XViewType, YViewType,
          false, -1, 1, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
        functor_type functor (alpha, A, x, beta, y);
        Kokkos::parallel_reduce (range, functor);
      }
      else { // beta != 0 && beta != 1
        typedef SingleLevelTransposeGEMV<AViewType, XViewType, YViewType,
          false, -1, -1, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
        functor_type functor (alpha, A, x, beta, y);
        Kokkos::parallel_reduce (range, functor);
      }
    }
  }
  else if (tr == 'C' || tr == 'c' || tr == 'H' || tr == 'h') { // conj xpose
    if (alpha == Kokkos::Details::ArithTraits<AlphaCoeffType>::zero ()) {
      if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::zero ()) {
        // Fill y with zeros
        Kokkos::deep_copy (y, Kokkos::Details::ArithTraits<y_value_type>::zero ());
      }
      else if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::one ()) {
        // Do nothing (y := 1 * y)
      }
      else { // beta != 0 && beta != 1
        typedef SingleLevelTransposeGEMV<AViewType, XViewType, YViewType,
          true, 0, -1, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
        functor_type functor (alpha, A, x, beta, y);
        Kokkos::parallel_reduce (range, functor);
      }
    }
    else if (alpha == Kokkos::Details::ArithTraits<AlphaCoeffType>::one ()) {
      if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::zero ()) {
        typedef SingleLevelTransposeGEMV<AViewType, XViewType, YViewType,
          true, 1, 0, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
        functor_type functor (alpha, A, x, beta, y);
        Kokkos::parallel_reduce (range, functor);
      }
      else if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::one ()) {
        typedef SingleLevelTransposeGEMV<AViewType, XViewType, YViewType,
          true, 1, 1, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
        functor_type functor (alpha, A, x, beta, y);
        Kokkos::parallel_reduce (range, functor);
      }
      else { // beta != 0 && beta != 1
        typedef SingleLevelTransposeGEMV<AViewType, XViewType, YViewType,
          true, 1, -1, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
        functor_type functor (alpha, A, x, beta, y);
        Kokkos::parallel_reduce (range, functor);
      }
    }
    else { // alpha != 0 and alpha != 1
      if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::zero ()) {
        typedef SingleLevelTransposeGEMV<AViewType, XViewType, YViewType,
          true, -1, 0, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
        functor_type functor (alpha, A, x, beta, y);
        Kokkos::parallel_reduce (range, functor);
      }
      else if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::one ()) {
        typedef SingleLevelTransposeGEMV<AViewType, XViewType, YViewType,
          true, -1, 1, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
        functor_type functor (alpha, A, x, beta, y);
        Kokkos::parallel_reduce (range, functor);
      }
      else { // beta != 0 && beta != 1
        typedef SingleLevelTransposeGEMV<AViewType, XViewType, YViewType,
          true, -1, -1, AlphaCoeffType, BetaCoeffType, IndexType> functor_type;
        functor_type functor (alpha, A, x, beta, y);
        Kokkos::parallel_reduce (range, functor);
      }
    }
  }
}

// Implementation of KokkosBlas::gemv.
template<class AViewType,
         class XViewType,
         class YViewType,
         class AlphaCoeffType,
         class BetaCoeffType>
struct GEMV {
  static void
  gemv (const char trans[],
        const AlphaCoeffType& alpha,
        const AViewType& A,
        const XViewType& x,
        const BetaCoeffType& beta,
        const YViewType& y)
  {
    static_assert (Kokkos::Impl::is_view<AViewType>::value,
                   "AViewType must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XViewType>::value,
                   "XViewType must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YViewType>::value,
                   "YViewType must be a Kokkos::View.");
    static_assert (static_cast<int> (AViewType::rank) == 2,
                   "AViewType must have rank 2.");
    static_assert (static_cast<int> (XViewType::rank) == 1,
                   "XViewType must have rank 1.");
    static_assert (static_cast<int> (YViewType::rank) == 1,
                   "YViewType must have rank 1.");

    typedef typename AViewType::size_type size_type;
    const size_type numRows = A.dimension_0 ();
    const size_type numCols = A.dimension_1 ();

    // Prefer int as the index type, but use a larger type if needed.
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numCols < static_cast<size_type> (INT_MAX)) {
      singleLevelGemv<AViewType, XViewType, YViewType, AlphaCoeffType,
        BetaCoeffType, int> (trans, alpha, A, x, beta, y);
    }
    else {
      singleLevelGemv<AViewType, XViewType, YViewType, AlphaCoeffType,
        BetaCoeffType, size_type> (trans, alpha, A, x, beta, y);
    }
  }
};

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::GEMV.  This is NOT for users!!!  All the
// declarations of full specializations go in this header file.  We
// may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//

#define KOKKOSBLAS_IMPL_GEMV_DECL( SCALAR, MATRIX_LAYOUT, VECTOR_LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template<> \
struct GEMV<Kokkos::View<const SCALAR**, \
                         MATRIX_LAYOUT, \
                         Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
            Kokkos::View<const SCALAR*, \
                         VECTOR_LAYOUT, \
                         Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
            Kokkos::View<SCALAR*, \
                         VECTOR_LAYOUT, \
                         Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
            SCALAR, \
            SCALAR> \
{ \
  typedef Kokkos::View<const SCALAR**, \
    MATRIX_LAYOUT, \
    Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > AVT; \
  typedef Kokkos::View<const SCALAR*, \
    VECTOR_LAYOUT, \
    Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > XVT; \
  typedef Kokkos::View<SCALAR*, \
    VECTOR_LAYOUT, \
    Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > YVT; \
  typedef SCALAR coeff_type; \
  \
  static void \
  gemv (const char trans[], \
        const coeff_type& alpha, \
        const AVT& A, \
        const XVT& x, \
        const coeff_type& beta, \
        const YVT& y); \
};

//
// Declarations of full specialization of KokkosBlas::Impl::GEMV.
// Their definitions go in .cpp file(s) in this source directory.
//

#ifdef KOKKOS_HAVE_SERIAL

KOKKOSBLAS_IMPL_GEMV_DECL( int, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace )
KOKKOSBLAS_IMPL_GEMV_DECL( long, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace )
KOKKOSBLAS_IMPL_GEMV_DECL( double, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace )

#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_OPENMP

// KOKKOSBLAS_IMPL_GEMV_DECL( int, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace )
// KOKKOSBLAS_IMPL_GEMV_DECL( long, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace )
// KOKKOSBLAS_IMPL_GEMV_DECL( double, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace )

#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD

KOKKOSBLAS_IMPL_GEMV_DECL( int, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace )
KOKKOSBLAS_IMPL_GEMV_DECL( long, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace )
KOKKOSBLAS_IMPL_GEMV_DECL( double, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace )

#endif // KOKKOS_HAVE_PTHREAD

#ifdef KOKKOS_HAVE_CUDA

KOKKOSBLAS_IMPL_GEMV_DECL( int, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace )
KOKKOSBLAS_IMPL_GEMV_DECL( long, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace )
KOKKOSBLAS_IMPL_GEMV_DECL( double, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace )

#endif // KOKKOS_HAVE_CUDA

//
// Macro for declarations of full specialization of
// KokkosBlas::Impl::GEMV.  This is NOT for users!!!
//

#define KOKKOSBLAS_IMPL_GEMV_DEF( SCALAR, MATRIX_LAYOUT, VECTOR_LAYOUT, EXEC_SPACE, MEM_SPACE ) \
void \
GEMV<Kokkos::View<const SCALAR**, \
                  MATRIX_LAYOUT,                                \
                  Kokkos::Device<EXEC_SPACE, MEM_SPACE>,            \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >,        \
     Kokkos::View<const SCALAR*,                                    \
                  VECTOR_LAYOUT,                                    \
                  Kokkos::Device<EXEC_SPACE, MEM_SPACE>,            \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >,        \
     Kokkos::View<SCALAR*,                                          \
                  VECTOR_LAYOUT,                                    \
                  Kokkos::Device<EXEC_SPACE, MEM_SPACE>,            \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >,        \
     SCALAR,                                                        \
     SCALAR>::                                                      \
gemv (const char trans[], \
      const coeff_type& alpha, \
      const AVT& A, \
      const XVT& x, \
      const coeff_type& beta, \
      const YVT& y) \
{ \
  static_assert (Kokkos::Impl::is_view<AVT>::value, \
                 "AVT must be a Kokkos::View."); \
  static_assert (Kokkos::Impl::is_view<XVT>::value, \
                 "XVT must be a Kokkos::View."); \
  static_assert (Kokkos::Impl::is_view<YVT>::value, \
                 "YVT must be a Kokkos::View."); \
  static_assert (static_cast<int> (AVT::rank) == 2, \
                 "AVT must have rank 2."); \
  static_assert (static_cast<int> (XVT::rank) == 1, \
                 "XVT must have rank 1."); \
  static_assert (static_cast<int> (YVT::rank) == 1, \
                 "YVT must have rank 1."); \
  \
  typedef AVT::size_type size_type; \
  typedef SCALAR coeff_type; \
  const size_type numRows = A.dimension_0 (); \
  const size_type numCols = A.dimension_1 (); \
  \
  if (numRows < static_cast<size_type> (INT_MAX) && \
      numCols < static_cast<size_type> (INT_MAX)) { \
    singleLevelGemv<AVT, XVT, YVT, coeff_type, \
      coeff_type, int> (trans, alpha, A, x, beta, y); \
  } \
  else { \
    singleLevelGemv<AVT, XVT, YVT, coeff_type, \
      coeff_type, size_type> (trans, alpha, A, x, beta, y); \
  } \
}

} // namespace Impl
} // namespace KokkosBlas

#endif // KOKKOS_BLAS2_MV_IMPL_GEMV_HPP_
