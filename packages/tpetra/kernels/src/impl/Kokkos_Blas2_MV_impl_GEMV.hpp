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

#ifdef HAVE_TPETRAKERNELS_MKL
#  include "mkl.h"
#endif // HAVE_TPETRAKERNELS_MKL
#ifdef KOKKOS_HAVE_OMP
#  include "omp.h"
#endif // KOKKOS_HAVE_OMP

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

#if ! defined(HAVE_TPETRAKERNELS_MKL)

template<class AViewType,
         class XViewType,
         class YViewType,
         class AlphaCoeffType,
         class BetaCoeffType,
         class IndexType,
         class ValueType = typename AViewType::non_const_value_type,
         const bool valueTypeSupported =
           (std::is_same<ValueType, double>::value ||
            std::is_same<ValueType, float>::value ||
            std::is_same<ValueType, ::Kokkos::complex<float> >::value ||
            std::is_same<ValueType, ::Kokkos::complex<double> >::value),
         const bool allValueTypesSame =
           (std::is_same<ValueType, typename XViewType::non_const_value_type>::value &&
            std::is_same<ValueType, typename YViewType::non_const_value_type>::value &&
            std::is_same<ValueType, AlphaCoeffType>::value &&
            std::is_same<ValueType, BetaCoeffType>::value),
         const bool layoutSupported =
           (std::is_same<typename AViewType::array_layout, ::Kokkos::LayoutLeft>::value ||
            std::is_same<typename AViewType::array_layout, ::Kokkos::LayoutRight>::value),
         const bool allLayoutsSame =
           (std::is_same<typename AViewType::array_layout, typename XViewType::array_layout>::value &&
            std::is_same<typename AViewType::array_layout, typename YViewType::array_layout>::value),
         const bool indexTypeSupported = std::is_same<IndexType, int>::value>
struct TryMklCblas
{
  static bool implemented () { return false; }

  static bool
  gemv (const char /* trans */ [],
        const AlphaCoeffType& /* alpha */,
        const AViewType& /* A */,
        const XViewType& /* x */,
        const BetaCoeffType& /* beta */,
        const YViewType& /* y */)
  {
    return false;
  }
};

#else

// Attempt to invoke MKL for gemv.  If it does actually call the MKL,
// return true, else return false.
template<class AViewType,
         class XViewType,
         class YViewType,
         class AlphaCoeffType,
         class BetaCoeffType,
         class IndexType,
         class ValueType = typename AViewType::non_const_value_type,
         const bool valueTypeSupported =
           (std::is_same<ValueType, double>::value ||
            std::is_same<ValueType, float>::value ||
            std::is_same<ValueType, ::Kokkos::complex<float> >::value ||
            std::is_same<ValueType, ::Kokkos::complex<double> >::value),
         const bool allValueTypesSame =
           (std::is_same<ValueType, typename XViewType::non_const_value_type>::value &&
            std::is_same<ValueType, typename YViewType::non_const_value_type>::value &&
            std::is_same<ValueType, AlphaCoeffType>::value &&
            std::is_same<ValueType, BetaCoeffType>::value),
         const bool layoutSupported =
           (std::is_same<typename AViewType::array_layout, ::Kokkos::LayoutLeft>::value ||
            std::is_same<typename AViewType::array_layout, ::Kokkos::LayoutRight>::value),
         const bool allLayoutsSame =
           (std::is_same<typename AViewType::array_layout, typename XViewType::array_layout>::value &&
            std::is_same<typename AViewType::array_layout, typename YViewType::array_layout>::value) >
struct TryMklCblas
{
  static bool implemented () { return false; }

  static bool
  gemv (const char /* trans */ [],
        const AlphaCoeffType& /* alpha */,
        const AViewType& /* A */,
        const XViewType& /* x */,
        const BetaCoeffType& /* beta */,
        const YViewType& /* y */)
  {
    return false;
  }
};

template<class AViewType,
         class XViewType,
         class YViewType,
         class AlphaCoeffType,
         class BetaCoeffType,
         class IndexType>
struct TryMklCblas<AViewType,
                   XViewType,
                   YViewType,
                   AlphaCoeffType,
                   BetaCoeffType,
                   IndexType,
                   double, true, true, true, true>
{
  static bool implemented () { return true; }

  static bool
  gemv (const char trans[],
        const AlphaCoeffType& alpha,
        const AViewType& A,
        const XViewType& x,
        const BetaCoeffType& beta,
        const YViewType& y)
  {
    typedef typename AViewType::array_layout array_layout;

    const MKL_INT numRows = static_cast<MKL_INT> (A.dimension_0 ());
    const MKL_INT numCols = static_cast<MKL_INT> (A.dimension_1 ());
    const CBLAS_LAYOUT layout =
      std::is_same<array_layout, Kokkos::LayoutRight>::value ?
      CblasRowMajor : CblasColMajor;
    const bool transpose = (trans[0] != 'N' && trans[0] != 'n');
    const bool conjugate = (trans[0] == 'C' || trans[0] == 'c' ||
                            trans[0] == 'H' || trans[0] == 'h');
    CBLAS_TRANSPOSE transCblas = CblasNoTrans;
    if (transpose) {
      if (conjugate) {
        transCblas = CblasConjTrans;
      }
      else {
        transCblas = CblasTrans;
      }
    }

    // FIXME (mfh 10 May 2016) I haven't quite figured out how to get
    // the strides right for the case of zero rows or columns, so I
    // handle these cases specially.
    if (numRows == 0) {
      if (transpose) {
        // Treat this like alpha == 0.
        if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::zero ()) {
          Kokkos::deep_copy (y, Kokkos::Details::ArithTraits<typename AViewType::non_const_value_type>::zero ());
        }
        else if (beta != Kokkos::Details::ArithTraits<BetaCoeffType>::one ()) {
          Kokkos::RangePolicy<typename AViewType::execution_space, IndexType>
            range (0, y.dimension_0 ());
          Kokkos::parallel_for (range, KOKKOS_LAMBDA (const IndexType& i) {
              y(i) = beta * y(i);
            });
        }
      }
      else {
        return true; // output vector y has length zero; nothing to do
      }
    }
    else if (numCols == 0) {
      if (transpose) {
        return true; // output vector y has length zero; nothing to do
      }
      else {
        // Treat this like alpha == 0.
        if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::zero ()) {
          Kokkos::deep_copy (y, Kokkos::Details::ArithTraits<typename AViewType::non_const_value_type>::zero ());
        }
        else if (beta != Kokkos::Details::ArithTraits<BetaCoeffType>::one ()) {
          Kokkos::RangePolicy<typename AViewType::execution_space, IndexType>
            range (0, y.dimension_0 ());
          Kokkos::parallel_for (range, KOKKOS_LAMBDA (const IndexType& i) {
              y(i) = beta * y(i);
            });
        }
      }
    }

    MKL_INT strides[8];
    A.stride (strides);
    MKL_INT lda;
    if (std::is_same<array_layout, Kokkos::LayoutLeft>::value) { // column major
      lda = (numCols > static_cast<MKL_INT> (1)) ? strides[1] : numRows;
      if (lda == 0) {
        lda = 1; // as the BLAS requires: lda >= max(1,m)
      }
    }
    else { // row major
      lda = (numRows > static_cast<MKL_INT> (1)) ? strides[0] : numCols;
      if (lda == 0) {
        lda = 1;
      }
    }
    const MKL_INT incx = 1;
    const MKL_INT incy = 1;

    cblas_dgemv (layout, transCblas, numRows, numCols,
                 alpha, A.ptr_on_device (), lda,
                 x.ptr_on_device (), incx,
                 beta, y.ptr_on_device (), incy);
    return true;
  }
};

template<class AViewType,
         class XViewType,
         class YViewType,
         class AlphaCoeffType,
         class BetaCoeffType,
         class IndexType>
struct TryMklCblas<AViewType,
                   XViewType,
                   YViewType,
                   AlphaCoeffType,
                   BetaCoeffType,
                   IndexType,
                   float, true, true, true, true>
{
  static bool implemented () { return true; }

  static bool
  gemv (const char trans[],
        const AlphaCoeffType& alpha,
        const AViewType& A,
        const XViewType& x,
        const BetaCoeffType& beta,
        const YViewType& y)
  {
    typedef typename AViewType::array_layout array_layout;

    const MKL_INT numRows = static_cast<MKL_INT> (A.dimension_0 ());
    const MKL_INT numCols = static_cast<MKL_INT> (A.dimension_1 ());
    const CBLAS_LAYOUT layout =
      std::is_same<array_layout, Kokkos::LayoutRight>::value ?
      CblasRowMajor : CblasColMajor;
    const bool transpose = (trans[0] != 'N' && trans[0] != 'n');
    const bool conjugate = (trans[0] == 'C' || trans[0] == 'c' ||
                            trans[0] == 'H' || trans[0] == 'h');
    CBLAS_TRANSPOSE transCblas = CblasNoTrans;
    if (transpose) {
      if (conjugate) {
        transCblas = CblasConjTrans;
      }
      else {
        transCblas = CblasTrans;
      }
    }

    // FIXME (mfh 10 May 2016) I haven't quite figured out how to get
    // the strides right for the case of zero rows or columns, so I
    // handle these cases specially.
    if (numRows == 0) {
      if (transpose) {
        // Treat this like alpha == 0.
        if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::zero ()) {
          Kokkos::deep_copy (y, Kokkos::Details::ArithTraits<typename AViewType::non_const_value_type>::zero ());
        }
        else if (beta != Kokkos::Details::ArithTraits<BetaCoeffType>::one ()) {
          Kokkos::RangePolicy<typename AViewType::execution_space, IndexType>
            range (0, y.dimension_0 ());
          Kokkos::parallel_for (range, KOKKOS_LAMBDA (const IndexType& i) {
              y(i) = beta * y(i);
            });
        }
      }
      else {
        return true; // output vector y has length zero; nothing to do
      }
    }
    else if (numCols == 0) {
      if (transpose) {
        return true; // output vector y has length zero; nothing to do
      }
      else {
        // Treat this like alpha == 0.
        if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::zero ()) {
          Kokkos::deep_copy (y, Kokkos::Details::ArithTraits<typename AViewType::non_const_value_type>::zero ());
        }
        else if (beta != Kokkos::Details::ArithTraits<BetaCoeffType>::one ()) {
          Kokkos::RangePolicy<typename AViewType::execution_space, IndexType>
            range (0, y.dimension_0 ());
          Kokkos::parallel_for (range, KOKKOS_LAMBDA (const IndexType& i) {
              y(i) = beta * y(i);
            });
        }
      }
    }

    MKL_INT strides[8];
    A.stride (strides);
    MKL_INT lda;
    if (std::is_same<array_layout, Kokkos::LayoutLeft>::value) { // column major
      lda = (numCols > static_cast<MKL_INT> (1)) ? strides[1] : numRows;
      if (lda == 0) {
        lda = 1; // as the BLAS requires: lda >= max(1,m)
      }
    }
    else { // row major
      lda = (numRows > static_cast<MKL_INT> (1)) ? strides[0] : numCols;
      if (lda == 0) {
        lda = 1;
      }
    }
    const MKL_INT incx = 1;
    const MKL_INT incy = 1;

    cblas_sgemv (layout, transCblas, numRows, numCols,
                 alpha, A.ptr_on_device (), lda,
                 x.ptr_on_device (), incx,
                 beta, y.ptr_on_device (), incy);
    return true;
  }
};

template<class AViewType,
         class XViewType,
         class YViewType,
         class AlphaCoeffType,
         class BetaCoeffType,
         class IndexType>
struct TryMklCblas<AViewType,
                   XViewType,
                   YViewType,
                   AlphaCoeffType,
                   BetaCoeffType,
                   IndexType,
                   ::Kokkos::complex<float>, true, true, true, true>
{
  static bool implemented () { return true; }

  static bool
  gemv (const char trans[],
        const AlphaCoeffType& alpha,
        const AViewType& A,
        const XViewType& x,
        const BetaCoeffType& beta,
        const YViewType& y)
  {
    typedef typename AViewType::array_layout array_layout;

    const CBLAS_LAYOUT layout =
      std::is_same<array_layout, Kokkos::LayoutRight>::value ?
      CblasRowMajor : CblasColMajor;

    const bool transpose = (trans[0] != 'N' && trans[0] != 'n');
    const bool conjugate = (trans[0] == 'C' || trans[0] == 'c' ||
                            trans[0] == 'H' || trans[0] == 'h');
    CBLAS_TRANSPOSE transCblas = CblasNoTrans;
    if (transpose) {
      if (conjugate) {
        transCblas = CblasConjTrans;
      }
      else {
        transCblas = CblasTrans;
      }
    }

    const MKL_INT numRows = static_cast<MKL_INT> (A.dimension_0 ());
    const MKL_INT numCols = static_cast<MKL_INT> (A.dimension_1 ());

    // FIXME (mfh 10 May 2016) I haven't quite figured out how to get
    // the strides right for the case of zero rows or columns, so I
    // handle these cases specially.
    if (numRows == 0) {
      if (transpose) {
        // Treat this like alpha == 0.
        if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::zero ()) {
          Kokkos::deep_copy (y, Kokkos::Details::ArithTraits<typename AViewType::non_const_value_type>::zero ());
        }
        else if (beta != Kokkos::Details::ArithTraits<BetaCoeffType>::one ()) {
          Kokkos::RangePolicy<typename AViewType::execution_space, IndexType>
            range (0, y.dimension_0 ());
          Kokkos::parallel_for (range, KOKKOS_LAMBDA (const IndexType& i) {
              y(i) = beta * y(i);
            });
        }
      }
      else {
        return true; // output vector y has length zero; nothing to do
      }
    }
    else if (numCols == 0) {
      if (transpose) {
        return true; // output vector y has length zero; nothing to do
      }
      else {
        // Treat this like alpha == 0.
        if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::zero ()) {
          Kokkos::deep_copy (y, Kokkos::Details::ArithTraits<typename AViewType::non_const_value_type>::zero ());
        }
        else if (beta != Kokkos::Details::ArithTraits<BetaCoeffType>::one ()) {
          Kokkos::RangePolicy<typename AViewType::execution_space, IndexType>
            range (0, y.dimension_0 ());
          Kokkos::parallel_for (range, KOKKOS_LAMBDA (const IndexType& i) {
              y(i) = beta * y(i);
            });
        }
      }
    }

    MKL_INT strides[8];
    A.stride (strides);
    MKL_INT lda;
    if (std::is_same<array_layout, Kokkos::LayoutLeft>::value) { // column major
      lda = (numCols > static_cast<MKL_INT> (1)) ? strides[1] : numRows;
      if (lda == 0) {
        lda = 1; // as the BLAS requires: lda >= max(1,m)
      }
    }
    else { // row major
      lda = (numRows > static_cast<MKL_INT> (1)) ? strides[0] : numCols;
      if (lda == 0) {
        lda = 1;
      }
    }
    const MKL_INT incx = 1;
    const MKL_INT incy = 1;

    // The complex versions of this routine take the coefficients by pointer.
    cblas_cgemv (layout, transCblas, numRows, numCols,
                 &alpha, A.ptr_on_device (), lda,
                 x.ptr_on_device (), incx,
                 &beta, y.ptr_on_device (), incy);
    return true;
  }
};

template<class AViewType,
         class XViewType,
         class YViewType,
         class AlphaCoeffType,
         class BetaCoeffType,
         class IndexType>
struct TryMklCblas<AViewType,
                   XViewType,
                   YViewType,
                   AlphaCoeffType,
                   BetaCoeffType,
                   IndexType,
                   ::Kokkos::complex<double>, true, true, true, true>
{
  static bool implemented () { return true; }

  static bool
  gemv (const char trans[],
        const AlphaCoeffType& alpha,
        const AViewType& A,
        const XViewType& x,
        const BetaCoeffType& beta,
        const YViewType& y)
  {
    typedef typename AViewType::array_layout array_layout;

    const CBLAS_LAYOUT layout =
      std::is_same<array_layout, Kokkos::LayoutRight>::value ?
      CblasRowMajor : CblasColMajor;

    const bool transpose = (trans[0] != 'N' && trans[0] != 'n');
    const bool conjugate = (trans[0] == 'C' || trans[0] == 'c' ||
                            trans[0] == 'H' || trans[0] == 'h');
    CBLAS_TRANSPOSE transCblas = CblasNoTrans;
    if (transpose) {
      if (conjugate) {
        transCblas = CblasConjTrans;
      }
      else {
        transCblas = CblasTrans;
      }
    }

    const MKL_INT numRows = static_cast<MKL_INT> (A.dimension_0 ());
    const MKL_INT numCols = static_cast<MKL_INT> (A.dimension_1 ());

    // FIXME (mfh 10 May 2016) I haven't quite figured out how to get
    // the strides right for the case of zero rows or columns, so I
    // handle these cases specially.
    if (numRows == 0) {
      if (transpose) {
        // Treat this like alpha == 0.
        if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::zero ()) {
          Kokkos::deep_copy (y, Kokkos::Details::ArithTraits<typename AViewType::non_const_value_type>::zero ());
        }
        else if (beta != Kokkos::Details::ArithTraits<BetaCoeffType>::one ()) {
          Kokkos::RangePolicy<typename AViewType::execution_space, IndexType>
            range (0, y.dimension_0 ());
          Kokkos::parallel_for (range, KOKKOS_LAMBDA (const IndexType& i) {
              y(i) = beta * y(i);
            });
        }
      }
      else {
        return true; // output vector y has length zero; nothing to do
      }
    }
    else if (numCols == 0) {
      if (transpose) {
        return true; // output vector y has length zero; nothing to do
      }
      else {
        // Treat this like alpha == 0.
        if (beta == Kokkos::Details::ArithTraits<BetaCoeffType>::zero ()) {
          Kokkos::deep_copy (y, Kokkos::Details::ArithTraits<typename AViewType::non_const_value_type>::zero ());
        }
        else if (beta != Kokkos::Details::ArithTraits<BetaCoeffType>::one ()) {
          Kokkos::RangePolicy<typename AViewType::execution_space, IndexType>
            range (0, y.dimension_0 ());
          Kokkos::parallel_for (range, KOKKOS_LAMBDA (const IndexType& i) {
              y(i) = beta * y(i);
            });
        }
      }
    }

    MKL_INT strides[8];
    A.stride (strides);
    MKL_INT lda;
    if (std::is_same<array_layout, Kokkos::LayoutLeft>::value) { // column major
      lda = (numCols > static_cast<MKL_INT> (1)) ? strides[1] : numRows;
      if (lda == 0) {
        lda = 1; // as the BLAS requires: lda >= max(1,m)
      }
    }
    else { // row major
      lda = (numRows > static_cast<MKL_INT> (1)) ? strides[0] : numCols;
      if (lda == 0) {
        lda = 1;
      }
    }
    const MKL_INT incx = 1;
    const MKL_INT incy = 1;

    // The complex versions of this routine take the coefficients by pointer.
    cblas_zgemv (layout, transCblas, numRows, numCols,
                 &alpha, A.ptr_on_device (), lda,
                 x.ptr_on_device (), incx,
                 &beta, y.ptr_on_device (), incy);
    return true;
  }
};

#endif // defined(HAVE_TPETRAKERNELS_MKL)

template<class AViewType,
         class XViewType,
         class YViewType,
         class AlphaCoeffType,
         class BetaCoeffType,
         class IndexType = typename AViewType::size_type>
bool
tryMklGemv (const char trans[],
            const AlphaCoeffType& alpha,
            const AViewType& A,
            const XViewType& x,
            const BetaCoeffType& beta,
            const YViewType& y)
{
  // FIXME (mfh 10 May 2016) See note below.  We use IndexType = int
  // explicitly here in order to avoid extra instantiations that don't
  // work anyway.  We should really use MKL_INT instead of int.
  typedef TryMklCblas<AViewType, XViewType, YViewType, AlphaCoeffType,
    BetaCoeffType, int> impl_type;
  typedef typename AViewType::size_type size_type;

  if (impl_type::implemented ()) {
    // FIXME (mfh 10 May 2016) MKL works for MKL_INT, which might not
    // necessarily be int.  However, assuming MKL_INT = int is safe
    // for now.
    if (A.dimension_0 () < static_cast<size_type> (INT_MAX) &&
        A.dimension_1 () < static_cast<size_type> (INT_MAX)) {
      return impl_type::gemv (trans, alpha, A, x, beta, y);
    }
    else {
      return false;
    }
  }
  else {
    return false;
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
