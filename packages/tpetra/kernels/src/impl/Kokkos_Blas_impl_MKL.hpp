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
#ifndef KOKKOS_BLAS_IMPL_MKL_HPP_
#define KOKKOS_BLAS_IMPL_MKL_HPP_

#include "TpetraKernels_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"

#ifdef HAVE_TPETRAKERNELS_MKL
#  include "mkl.h"
#endif // HAVE_TPETRAKERNELS_MKL

namespace KokkosBlas {
namespace Impl {

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

} // namespace Impl
} // namespace KokkosBlas

#endif // KOKKOS_BLAS_IMPL_MKL_HPP_
