//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

/// \file
/// \brief Interfaces for the Kokkos sparse-matrix-vector multiply
///

#ifndef KOKKOSSPARSE_SPMV_HPP_
#define KOKKOSSPARSE_SPMV_HPP_

#include "KokkosKernels_helpers.hpp"
#include "KokkosKernels_Controls.hpp"
#include "KokkosSparse_spmv_spec.hpp"
#include "KokkosSparse_spmv_struct_spec.hpp"
#include "KokkosSparse_spmv_bsrmatrix_spec.hpp"
#include <type_traits>
#include "KokkosSparse_BsrMatrix.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosBlas1_scal.hpp"
#include "KokkosKernels_Utils.hpp"
#include "KokkosKernels_Error.hpp"

namespace KokkosSparse {

namespace {
struct RANK_ONE {};
struct RANK_TWO {};
}  // namespace

/// \brief Tag-dispatch for \c Kokkos sparse matrix-vector multiply on single
/// vector
///
///
/// \tparam AMatrix  A KokkosSparse::CrsMatrix, or KokkosSparse::BsrMatrix
///
/// \param controls [in] kokkos-kernels control structure.
/// \param mode [in]
/// \param alpha [in] Scalar multiplier for the matrix A.
/// \param A [in] The sparse matrix A.
/// \param x [in] A vector.
/// \param beta [in] Scalar multiplier for the multivector y.
/// \param y [in/out] vector.
/// \param RANK_ONE tag dispatch
///
#ifdef DOXY  // documentation version
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
#else
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector,
          typename std::enable_if<
              KokkosSparse::is_crs_matrix<AMatrix>::value>::type* = nullptr>
#endif
void spmv(KokkosKernels::Experimental::Controls controls, const char mode[],
          const AlphaType& alpha, const AMatrix& A, const XVector& x,
          const BetaType& beta, const YVector& y, const RANK_ONE) {

  // Make sure that x and y have the same rank.
  static_assert(
      static_cast<int>(XVector::rank) == static_cast<int>(YVector::rank),
      "KokkosSparse::spmv: Vector ranks do not match.");
  // Make sure that x (and therefore y) is rank 1.
  static_assert(static_cast<int>(XVector::rank) == 1,
                "KokkosSparse::spmv: Both Vector inputs must have rank 1 "
                "in order to call this specialization of spmv.");
  // Make sure that y is non-const.
  static_assert(std::is_same<typename YVector::value_type,
                             typename YVector::non_const_value_type>::value,
                "KokkosSparse::spmv: Output Vector must be non-const.");

  // Check compatibility of dimensions at run time.
  if ((mode[0] == NoTranspose[0]) || (mode[0] == Conjugate[0])) {
    if ((x.extent(1) != y.extent(1)) ||
        (static_cast<size_t>(A.numCols()) > static_cast<size_t>(x.extent(0))) ||
        (static_cast<size_t>(A.numRows()) > static_cast<size_t>(y.extent(0)))) {
      std::ostringstream os;
      os << "KokkosSparse::spmv: Dimensions do not match: "
         << ", A: " << A.numRows() << " x " << A.numCols()
         << ", x: " << x.extent(0) << " x " << x.extent(1)
         << ", y: " << y.extent(0) << " x " << y.extent(1);
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  } else {
    if ((x.extent(1) != y.extent(1)) ||
        (static_cast<size_t>(A.numCols()) > static_cast<size_t>(y.extent(0))) ||
        (static_cast<size_t>(A.numRows()) > static_cast<size_t>(x.extent(0)))) {
      std::ostringstream os;
      os << "KokkosSparse::spmv: Dimensions do not match (transpose): "
         << ", A: " << A.numRows() << " x " << A.numCols()
         << ", x: " << x.extent(0) << " x " << x.extent(1)
         << ", y: " << y.extent(0) << " x " << y.extent(1);
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  }

  typedef KokkosSparse::CrsMatrix<
      typename AMatrix::const_value_type, typename AMatrix::const_ordinal_type,
      typename AMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>,
      typename AMatrix::const_size_type>
      AMatrix_Internal;

  typedef Kokkos::View<
      typename XVector::const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<XVector>::array_layout,
      typename XVector::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      XVector_Internal;

  typedef Kokkos::View<
      typename YVector::non_const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<YVector>::array_layout,
      typename YVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      YVector_Internal;

  AMatrix_Internal A_i = A;
  XVector_Internal x_i = x;
  YVector_Internal y_i = y;

  if (alpha == Kokkos::ArithTraits<AlphaType>::zero() || A_i.numRows() == 0 ||
      A_i.numCols() == 0 || A_i.nnz() == 0) {
    // This is required to maintain semantics of KokkosKernels native SpMV:
    // if y contains NaN but beta = 0, the result y should be filled with 0.
    // For example, this is useful for passing in uninitialized y and beta=0.
    if (beta == Kokkos::ArithTraits<BetaType>::zero())
      Kokkos::deep_copy(y_i, Kokkos::ArithTraits<BetaType>::zero());
    else
      KokkosBlas::scal(y_i, beta, y_i);
    return;
  }

  // Whether to call KokkosKernel's native implementation, even if a TPL impl is
  // available
  bool useFallback = controls.isParameter("algorithm") &&
                     controls.getParameter("algorithm") == "native";

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  // cuSPARSE does not support the conjugate mode (C), and cuSPARSE 9 only
  // supports the normal (N) mode.
  if (std::is_same<typename AMatrix_Internal::memory_space,
                   Kokkos::CudaSpace>::value ||
      std::is_same<typename AMatrix_Internal::memory_space,
                   Kokkos::CudaUVMSpace>::value) {
#if (9000 <= CUDA_VERSION)
    useFallback = useFallback || (mode[0] != NoTranspose[0]);
#endif
#if defined(CUSPARSE_VERSION) && (10300 <= CUSPARSE_VERSION)
    useFallback = useFallback || (mode[0] == Conjugate[0]);
#endif
  }
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
  if (std::is_same<typename AMatrix_Internal::memory_space,
                   Kokkos::Experimental::HIPSpace>::value) {
    useFallback = useFallback || (mode[0] != NoTranspose[0]);
  }
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
  if (std::is_same<typename AMatrix_Internal::memory_space,
                   Kokkos::HostSpace>::value) {
    useFallback = useFallback || (mode[0] == Conjugate[0]);
  }
#endif

  if (useFallback) {
    // Explicitly call the non-TPL SPMV implementation
    std::string label =
        "KokkosSparse::spmv[NATIVE," +
        Kokkos::ArithTraits<
            typename AMatrix_Internal::non_const_value_type>::name() +
        "]";
    Kokkos::Profiling::pushRegion(label);
    Impl::SPMV<typename AMatrix_Internal::value_type,
               typename AMatrix_Internal::ordinal_type,
               typename AMatrix_Internal::device_type,
               typename AMatrix_Internal::memory_traits,
               typename AMatrix_Internal::size_type,
               typename XVector_Internal::value_type*,
               typename XVector_Internal::array_layout,
               typename XVector_Internal::device_type,
               typename XVector_Internal::memory_traits,
               typename YVector_Internal::value_type*,
               typename YVector_Internal::array_layout,
               typename YVector_Internal::device_type,
               typename YVector_Internal::memory_traits, false>::spmv(controls,
                                                                      mode,
                                                                      alpha,
                                                                      A_i, x_i,
                                                                      beta,
                                                                      y_i);
    Kokkos::Profiling::popRegion();
  } else {
    // note: the cuSPARSE spmv wrapper defines a profiling region, so one is not
    // needed here.
    Impl::SPMV<typename AMatrix_Internal::value_type,
               typename AMatrix_Internal::ordinal_type,
               typename AMatrix_Internal::device_type,
               typename AMatrix_Internal::memory_traits,
               typename AMatrix_Internal::size_type,
               typename XVector_Internal::value_type*,
               typename XVector_Internal::array_layout,
               typename XVector_Internal::device_type,
               typename XVector_Internal::memory_traits,
               typename YVector_Internal::value_type*,
               typename YVector_Internal::array_layout,
               typename YVector_Internal::device_type,
               typename YVector_Internal::memory_traits>::spmv(controls, mode,
                                                               alpha, A_i, x_i,
                                                               beta, y_i);
  }
}

#ifdef DOXY  // hide SFINAE
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
#else
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector,
          typename std::enable_if<KokkosSparse::Experimental::is_bsr_matrix<
              AMatrix>::value>::type* = nullptr>
#endif
void spmv(KokkosKernels::Experimental::Controls controls, const char mode[],
          const AlphaType& alpha, const AMatrix& A, const XVector& x,
          const BetaType& beta, const YVector& y, const RANK_ONE) {
  // Make sure that x and y have the same rank.
  static_assert(
      static_cast<int>(XVector::rank) == static_cast<int>(YVector::rank),
      "KokkosSparse::spmv: Vector ranks do not match.");
  // Make sure that x (and therefore y) is rank 1.
  static_assert(static_cast<int>(XVector::rank) == 1,
                "KokkosSparse::spmv: Both Vector inputs must have rank 1 "
                "in order to call this specialization of spmv.");
  // Make sure that y is non-const.
  static_assert(std::is_same<typename YVector::value_type,
                             typename YVector::non_const_value_type>::value,
                "KokkosSparse::spmv: Output Vector must be non-const.");

  //
  if (A.blockDim() == 1) {
    KokkosSparse::CrsMatrix<
        typename AMatrix::value_type, typename AMatrix::ordinal_type,
        typename AMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>,
        typename AMatrix::size_type>
        Acrs("bsr_to_crs", A.numCols(), A.values, A.graph);
    KokkosSparse::spmv(controls, mode, alpha, Acrs, x, beta, y, RANK_ONE());
    return;
  }
  // Check compatibility of dimensions at run time.
  if ((mode[0] == NoTranspose[0]) || (mode[0] == Conjugate[0])) {
    if ((x.extent(1) != y.extent(1)) ||
        (static_cast<size_t>(A.numCols() * A.blockDim()) !=
         static_cast<size_t>(x.extent(0))) ||
        (static_cast<size_t>(A.numRows() * A.blockDim()) !=
         static_cast<size_t>(y.extent(0)))) {
      std::ostringstream os;
      os << "KokkosSparse::spmv (BsrMatrix): Dimensions do not match: "
         << ", A: " << A.numRows() * A.blockDim() << " x "
         << A.numCols() * A.blockDim() << ", x: " << x.extent(0) << " x "
         << x.extent(1) << ", y: " << y.extent(0) << " x " << y.extent(1);

      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  } else {
    if ((x.extent(1) != y.extent(1)) ||
        (static_cast<size_t>(A.numCols() * A.blockDim()) !=
         static_cast<size_t>(y.extent(0))) ||
        (static_cast<size_t>(A.numRows() * A.blockDim()) !=
         static_cast<size_t>(x.extent(0)))) {
      std::ostringstream os;
      os << "KokkosSparse::spmv (BsrMatrix): Dimensions do not match "
            "(transpose): "
         << ", A: " << A.numRows() * A.blockDim() << " x "
         << A.numCols() * A.blockDim() << ", x: " << x.extent(0) << " x "
         << x.extent(1) << ", y: " << y.extent(0) << " x " << y.extent(1);

      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  }
  //
  typedef KokkosSparse::Experimental::BsrMatrix<
      typename AMatrix::const_value_type, typename AMatrix::const_ordinal_type,
      typename AMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>,
      typename AMatrix::const_size_type>
      AMatrix_Internal;

  typedef Kokkos::View<
      typename XVector::const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<XVector>::array_layout,
      typename XVector::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      XVector_Internal;

  typedef Kokkos::View<
      typename YVector::non_const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<YVector>::array_layout,
      typename YVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      YVector_Internal;

  AMatrix_Internal A_i(A);
  XVector_Internal x_i(x);
  YVector_Internal y_i(y);

  if (alpha == Kokkos::ArithTraits<AlphaType>::zero() || A_i.numRows() == 0 ||
      A_i.numCols() == 0 || A_i.nnz() == 0) {
    // This is required to maintain semantics of KokkosKernels native SpMV:
    // if y contains NaN but beta = 0, the result y should be filled with 0.
    // For example, this is useful for passing in uninitialized y and beta=0.
    if (beta == Kokkos::ArithTraits<BetaType>::zero())
      Kokkos::deep_copy(y_i, Kokkos::ArithTraits<BetaType>::zero());
    else
      KokkosBlas::scal(y_i, beta, y_i);
    return;
  }

  //
  // Whether to call KokkosKernel's native implementation, even if a TPL impl is
  // available
  bool useFallback = controls.isParameter("algorithm") &&
                     controls.getParameter("algorithm") == "native";

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  // cuSPARSE does not support the modes (C), (T), (H)
  if (std::is_same<typename AMatrix_Internal::memory_space,
                   Kokkos::CudaSpace>::value ||
      std::is_same<typename AMatrix_Internal::memory_space,
                   Kokkos::CudaUVMSpace>::value) {
    useFallback = useFallback || (mode[0] != NoTranspose[0]);
  }
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
  if (std::is_same<typename AMatrix_Internal::memory_space,
                   Kokkos::HostSpace>::value) {
    useFallback = useFallback || (mode[0] == Conjugate[0]);
  }
#endif

  if (useFallback) {
    // Explicitly call the non-TPL SPMV_BSRMATRIX implementation
    std::string label =
        "KokkosSparse::spmv[NATIVE,BSRMATRIX," +
        Kokkos::ArithTraits<
            typename AMatrix_Internal::non_const_value_type>::name() +
        "]";
    Kokkos::Profiling::pushRegion(label);
    Experimental::Impl::SPMV_BSRMATRIX<
        typename AMatrix_Internal::const_value_type,
        typename AMatrix_Internal::const_ordinal_type,
        typename AMatrix_Internal::device_type,
        typename AMatrix_Internal::memory_traits,
        typename AMatrix_Internal::const_size_type,
        typename XVector_Internal::const_value_type*,
        typename XVector_Internal::array_layout,
        typename XVector_Internal::device_type,
        typename XVector_Internal::memory_traits,
        typename YVector_Internal::value_type*,
        typename YVector_Internal::array_layout,
        typename YVector_Internal::device_type,
        typename YVector_Internal::memory_traits,
        false>::spmv_bsrmatrix(controls, mode, alpha, A_i, x_i, beta, y_i);
    Kokkos::Profiling::popRegion();
  } else {
#define __SPMV_TYPES__                               \
  typename AMatrix_Internal::const_value_type,       \
      typename AMatrix_Internal::const_ordinal_type, \
      typename AMatrix_Internal::device_type,        \
      typename AMatrix_Internal::memory_traits,      \
      typename AMatrix_Internal::const_size_type,    \
      typename XVector_Internal::const_value_type*,  \
      typename XVector_Internal::array_layout,       \
      typename XVector_Internal::device_type,        \
      typename XVector_Internal::memory_traits,      \
      typename YVector_Internal::value_type*,        \
      typename YVector_Internal::array_layout,       \
      typename YVector_Internal::device_type,        \
      typename YVector_Internal::memory_traits

    constexpr bool tpl_spec_avail =
        KokkosSparse::Experimental::Impl::spmv_bsrmatrix_tpl_spec_avail<
            __SPMV_TYPES__>::value;

    constexpr bool eti_spec_avail =
        tpl_spec_avail
            ? KOKKOSKERNELS_IMPL_COMPILE_LIBRARY /* force FALSE in app/test */
            : KokkosSparse::Experimental::Impl::spmv_bsrmatrix_eti_spec_avail<
                  __SPMV_TYPES__>::value;

    Experimental::Impl::SPMV_BSRMATRIX<__SPMV_TYPES__, tpl_spec_avail,
                                       eti_spec_avail>::spmv_bsrmatrix(controls,
                                                                       mode,
                                                                       alpha,
                                                                       A_i, x_i,
                                                                       beta,
                                                                       y_i);

#undef __SPMV_TYPES__
  }
}

template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector, class XLayout = typename XVector::array_layout>
struct SPMV2D1D {
  static bool spmv2d1d(const char mode[], const AlphaType& alpha,
                       const AMatrix& A, const XVector& x, const BetaType& beta,
                       const YVector& y);
};

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || !defined(KOKKOSKERNELS_ETI_ONLY)
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
struct SPMV2D1D<AlphaType, AMatrix, XVector, BetaType, YVector,
                Kokkos::LayoutStride> {
  static bool spmv2d1d(const char mode[], const AlphaType& alpha,
                       const AMatrix& A, const XVector& x, const BetaType& beta,
                       const YVector& y) {
    spmv(mode, alpha, A, x, beta, y);
    return true;
  }
#else
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
struct SPMV2D1D<AlphaType, AMatrix, XVector, BetaType, YVector,
                Kokkos::LayoutStride> {
  static bool spmv2d1d(const char /*mode*/[], const AlphaType& /*alpha*/,
                       const AMatrix& /*A*/, const XVector& /*x*/,
                       const BetaType& /*beta*/, const YVector& /*y*/) {
    return false;
  }
#endif
};

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || !defined(KOKKOSKERNELS_ETI_ONLY)
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
struct SPMV2D1D<AlphaType, AMatrix, XVector, BetaType, YVector,
                Kokkos::LayoutLeft> {
  static bool spmv2d1d(const char mode[], const AlphaType& alpha,
                       const AMatrix& A, const XVector& x, const BetaType& beta,
                       const YVector& y) {
    spmv(mode, alpha, A, x, beta, y);
    return true;
  }
#else
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
struct SPMV2D1D<AlphaType, AMatrix, XVector, BetaType, YVector,
                Kokkos::LayoutLeft> {
  static bool spmv2d1d(const char /*mode*/[], const AlphaType& /*alpha*/,
                       const AMatrix& /*A*/, const XVector& /*x*/,
                       const BetaType& /*beta*/, const YVector& /*y*/) {
    return false;
  }
#endif
};

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || !defined(KOKKOSKERNELS_ETI_ONLY)
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
struct SPMV2D1D<AlphaType, AMatrix, XVector, BetaType, YVector,
                Kokkos::LayoutRight> {
  static bool spmv2d1d(const char mode[], const AlphaType& alpha,
                       const AMatrix& A, const XVector& x, const BetaType& beta,
                       const YVector& y) {
    spmv(mode, alpha, A, x, beta, y);
    return true;
  }
#else
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
struct SPMV2D1D<AlphaType, AMatrix, XVector, BetaType, YVector,
                Kokkos::LayoutRight> {
  static bool spmv2d1d(const char /*mode*/[], const AlphaType& /*alpha*/,
                       const AMatrix& /*A*/, const XVector& /*x*/,
                       const BetaType& /*beta*/, const YVector& /*y*/) {
    return false;
  }
#endif
};

/// \brief Tag-dispatch sparse matrix-vector multiply on multivectors
///
/// \tparam AMatrix A KokkosSparse::CrsMatrix,
/// KokkosSparse::Experimental::BsrMatrix
///
/// \param controls [in] kokkos-kernels control structure.
/// \param mode [in] \c "N" for no transpose
/// \param alpha [in] Scalar multiplier for the matrix A.
/// \param A [in] The sparse matrix A.
/// \param x [in] A multivector (rank-2 Kokkos::View).
/// \param beta [in] Scalar multiplier for the multivector y.
/// \param y [in/out] multivector (rank-2 Kokkos::View).
/// \param RANK_TWO tag-dispatch
///
#ifdef DOXY
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
#else
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector,
          typename std::enable_if<
              KokkosSparse::is_crs_matrix<AMatrix>::value>::type* = nullptr>
#endif
void spmv(KokkosKernels::Experimental::Controls controls, const char mode[],
          const AlphaType& alpha, const AMatrix& A, const XVector& x,
          const BetaType& beta, const YVector& y, const RANK_TWO) {

  // Make sure that x and y have the same rank.
  static_assert(
      static_cast<int>(XVector::rank) == static_cast<int>(YVector::rank),
      "KokkosSparse::spmv: Vector ranks do not match.");
  // Make sure that x (and therefore y) is rank 2.
  static_assert(static_cast<int>(XVector::rank) == 2,
                "KokkosSparse::spmv: Both Vector inputs must have rank 2 "
                "in order to call this specialization of spmv.");
  // Make sure that y is non-const.
  static_assert(std::is_same<typename YVector::value_type,
                             typename YVector::non_const_value_type>::value,
                "KokkosSparse::spmv: Output Vector must be non-const.");

  // Check compatibility of dimensions at run time.
  if ((mode[0] == NoTranspose[0]) || (mode[0] == Conjugate[0])) {
    if ((x.extent(1) != y.extent(1)) ||
        (static_cast<size_t>(A.numCols()) > static_cast<size_t>(x.extent(0))) ||
        (static_cast<size_t>(A.numRows()) > static_cast<size_t>(y.extent(0)))) {
      std::ostringstream os;
      os << "KokkosBlas::spmv: Dimensions do not match: "
         << ", A: " << A.numRows() << " x " << A.numCols()
         << ", x: " << x.extent(0) << " x " << x.extent(1)
         << ", y: " << y.extent(0) << " x " << y.extent(1);
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  } else {
    if ((x.extent(1) != y.extent(1)) ||
        (static_cast<size_t>(A.numCols()) > static_cast<size_t>(y.extent(0))) ||
        (static_cast<size_t>(A.numRows()) > static_cast<size_t>(x.extent(0)))) {
      std::ostringstream os;
      os << "KokkosBlas::spmv: Dimensions do not match (transpose): "
         << ", A: " << A.numRows() << " x " << A.numCols()
         << ", x: " << x.extent(0) << " x " << x.extent(1)
         << ", y: " << y.extent(0) << " x " << y.extent(1);
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  }

  typedef KokkosSparse::CrsMatrix<
      typename AMatrix::const_value_type, typename AMatrix::const_ordinal_type,
      typename AMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>,
      typename AMatrix::const_size_type>
      AMatrix_Internal;

  AMatrix_Internal A_i = A;

  // Call single-vector version if appropriate
  if (x.extent(1) == 1) {
    typedef Kokkos::View<
        typename XVector::const_value_type*,
        typename KokkosKernels::Impl::GetUnifiedLayout<XVector>::array_layout,
        typename XVector::device_type,
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
        XVector_SubInternal;
    typedef Kokkos::View<
        typename YVector::non_const_value_type*,
        typename KokkosKernels::Impl::GetUnifiedLayout<YVector>::array_layout,
        typename YVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
        YVector_SubInternal;

    XVector_SubInternal x_i = Kokkos::subview(x, Kokkos::ALL(), 0);
    YVector_SubInternal y_i = Kokkos::subview(y, Kokkos::ALL(), 0);

    // spmv (mode, alpha, A, x_i, beta, y_i);
    using impl_type = SPMV2D1D<AlphaType, AMatrix_Internal, XVector_SubInternal,
                               BetaType, YVector_SubInternal,
                               typename XVector_SubInternal::array_layout>;
    if (impl_type::spmv2d1d(mode, alpha, A, x_i, beta, y_i)) {
      return;
    }
  }
  {
    typedef Kokkos::View<
        typename XVector::const_value_type**, typename XVector::array_layout,
        typename XVector::device_type,
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
        XVector_Internal;

    typedef Kokkos::View<typename YVector::non_const_value_type**,
                         typename YVector::array_layout,
                         typename YVector::device_type,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >
        YVector_Internal;

    XVector_Internal x_i = x;
    YVector_Internal y_i = y;

    bool useNative = false;

// cusparseSpMM does not support conjugate mode
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
    useNative = useNative || (Conjugate[0] == mode[0]);
#endif
    useNative = useNative || (controls.isParameter("algorithm") &&
                              (controls.getParameter("algorithm") == "native"));

    if (useNative) {
      return Impl::SPMV_MV<
          typename AMatrix_Internal::value_type,
          typename AMatrix_Internal::ordinal_type,
          typename AMatrix_Internal::device_type,
          typename AMatrix_Internal::memory_traits,
          typename AMatrix_Internal::size_type,
          typename XVector_Internal::value_type**,
          typename XVector_Internal::array_layout,
          typename XVector_Internal::device_type,
          typename XVector_Internal::memory_traits,
          typename YVector_Internal::value_type**,
          typename YVector_Internal::array_layout,
          typename YVector_Internal::device_type,
          typename YVector_Internal::memory_traits,
          std::is_integral<typename AMatrix_Internal::value_type>::value,
          false>::spmv_mv(controls, mode, alpha, A_i, x_i, beta, y_i);
    } else {
      return Impl::SPMV_MV<
          typename AMatrix_Internal::value_type,
          typename AMatrix_Internal::ordinal_type,
          typename AMatrix_Internal::device_type,
          typename AMatrix_Internal::memory_traits,
          typename AMatrix_Internal::size_type,
          typename XVector_Internal::value_type**,
          typename XVector_Internal::array_layout,
          typename XVector_Internal::device_type,
          typename XVector_Internal::memory_traits,
          typename YVector_Internal::value_type**,
          typename YVector_Internal::array_layout,
          typename YVector_Internal::device_type,
          typename YVector_Internal::memory_traits>::spmv_mv(controls, mode,
                                                             alpha, A_i, x_i,
                                                             beta, y_i);
    }
  }
}

#ifdef DOXY  // hide SFINAE
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
#else
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector,
          typename std::enable_if<KokkosSparse::Experimental::is_bsr_matrix<
              AMatrix>::value>::type* = nullptr>
#endif
void spmv(KokkosKernels::Experimental::Controls controls, const char mode[],
          const AlphaType& alpha, const AMatrix& A, const XVector& x,
          const BetaType& beta, const YVector& y, const RANK_TWO) {
  // Make sure that x and y have the same rank.
  static_assert(
      static_cast<int>(XVector::rank) == static_cast<int>(YVector::rank),
      "KokkosSparse::spmv: Vector ranks do not match.");
  // Make sure that x (and therefore y) is rank 2.
  static_assert(static_cast<int>(XVector::rank) == 2,
                "KokkosSparse::spmv: Both Vector inputs must have rank 2 "
                "in order to call this specialization of spmv.");
  // Make sure that y is non-const.
  static_assert(std::is_same<typename YVector::value_type,
                             typename YVector::non_const_value_type>::value,
                "KokkosSparse::spmv: Output Vector must be non-const.");

  //
  if (A.blockDim() == 1) {
    KokkosSparse::CrsMatrix<
        typename AMatrix::value_type, typename AMatrix::ordinal_type,
        typename AMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>,
        typename AMatrix::size_type>
        Acrs("bsr_to_crs", A.numCols(), A.values, A.graph);
    KokkosSparse::spmv(controls, mode, alpha, Acrs, x, beta, y, RANK_TWO());
    return;
  }
  // Check compatibility of dimensions at run time.
  if ((mode[0] == NoTranspose[0]) || (mode[0] == Conjugate[0])) {
    if ((x.extent(1) != y.extent(1)) ||
        (static_cast<size_t>(A.numCols() * A.blockDim()) !=
         static_cast<size_t>(x.extent(0))) ||
        (static_cast<size_t>(A.numRows() * A.blockDim()) !=
         static_cast<size_t>(y.extent(0)))) {
      std::ostringstream os;
      os << "KokkosSparse::spmv (BsrMatrix): Dimensions do not match: "
         << ", A: " << A.numRows() * A.blockDim() << " x "
         << A.numCols() * A.blockDim() << ", x: " << x.extent(0) << " x "
         << x.extent(1) << ", y: " << y.extent(0) << " x " << y.extent(1);

      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  } else {
    if ((x.extent(1) != y.extent(1)) ||
        (static_cast<size_t>(A.numCols() * A.blockDim()) !=
         static_cast<size_t>(y.extent(0))) ||
        (static_cast<size_t>(A.numRows() * A.blockDim()) !=
         static_cast<size_t>(x.extent(0)))) {
      std::ostringstream os;
      os << "KokkosSparse::spmv (BsrMatrix): Dimensions do not match "
            "(transpose): "
         << ", A: " << A.numRows() * A.blockDim() << " x "
         << A.numCols() * A.blockDim() << ", x: " << x.extent(0) << " x "
         << x.extent(1) << ", y: " << y.extent(0) << " x " << y.extent(1);

      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  }
  //
  typedef KokkosSparse::Experimental::BsrMatrix<
      typename AMatrix::const_value_type, typename AMatrix::const_ordinal_type,
      typename AMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>,
      typename AMatrix::const_size_type>
      AMatrix_Internal;
  AMatrix_Internal A_i(A);

  typedef Kokkos::View<
      typename XVector::const_value_type**,
      typename KokkosKernels::Impl::GetUnifiedLayout<XVector>::array_layout,
      typename XVector::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      XVector_Internal;
  XVector_Internal x_i(x);

  typedef Kokkos::View<
      typename YVector::non_const_value_type**,
      typename KokkosKernels::Impl::GetUnifiedLayout<YVector>::array_layout,
      typename YVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      YVector_Internal;
  YVector_Internal y_i(y);
  //
  if (alpha == Kokkos::ArithTraits<AlphaType>::zero() || A_i.numRows() == 0 ||
      A_i.numCols() == 0 || A_i.nnz() == 0) {
    // This is required to maintain semantics of KokkosKernels native SpMV:
    // if y contains NaN but beta = 0, the result y should be filled with 0.
    // For example, this is useful for passing in uninitialized y and beta=0.
    if (beta == Kokkos::ArithTraits<BetaType>::zero())
      Kokkos::deep_copy(y_i, Kokkos::ArithTraits<BetaType>::zero());
    else
      KokkosBlas::scal(y_i, beta, y_i);
    return;
  }
  //
  // Call single-vector version if appropriate
  //
  if (x.extent(1) == 1) {
    typedef Kokkos::View<
        typename XVector::const_value_type*,
        typename KokkosKernels::Impl::GetUnifiedLayout<XVector>::array_layout,
        typename XVector::device_type,
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
        XVector_SubInternal;
    typedef Kokkos::View<
        typename YVector::non_const_value_type*,
        typename KokkosKernels::Impl::GetUnifiedLayout<YVector>::array_layout,
        typename YVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
        YVector_SubInternal;

    XVector_SubInternal x_0 = Kokkos::subview(x_i, Kokkos::ALL(), 0);
    YVector_SubInternal y_0 = Kokkos::subview(y_i, Kokkos::ALL(), 0);

    return spmv(controls, mode, alpha, A_i, x_0, beta, y_0, RANK_ONE());
  }
  //
  // Whether to call KokkosKernel's native implementation, even if a TPL impl is
  // available
  bool useFallback =
      controls.isParameter("algorithm") &&
      (controls.getParameter("algorithm") == "native" ||
       controls.getParameter("algorithm") == "experimental_bsr_tc");

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  // cuSPARSE does not support the modes (C), (T), (H)
  if (std::is_same<typename AMatrix_Internal::memory_space,
                   Kokkos::CudaSpace>::value ||
      std::is_same<typename AMatrix_Internal::memory_space,
                   Kokkos::CudaUVMSpace>::value) {
    useFallback = useFallback || (mode[0] != NoTranspose[0]);
  }
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
  if (std::is_same<typename AMatrix_Internal::memory_space,
                   Kokkos::HostSpace>::value) {
    useFallback = useFallback || (mode[0] == Conjugate[0]);
  }
#endif

  if (useFallback) {
    // Explicitly call the non-TPL SPMV_BSRMATRIX implementation
    std::string label =
        "KokkosSparse::spmv[NATIVE,BSMATRIX," +
        Kokkos::ArithTraits<
            typename AMatrix_Internal::non_const_value_type>::name() +
        "]";
    Kokkos::Profiling::pushRegion(label);
    Experimental::Impl::SPMV_MV_BSRMATRIX<
        typename AMatrix_Internal::const_value_type,
        typename AMatrix_Internal::const_ordinal_type,
        typename AMatrix_Internal::device_type,
        typename AMatrix_Internal::memory_traits,
        typename AMatrix_Internal::const_size_type,
        typename XVector_Internal::const_value_type**,
        typename XVector_Internal::array_layout,
        typename XVector_Internal::device_type,
        typename XVector_Internal::memory_traits,
        typename YVector_Internal::value_type**,
        typename YVector_Internal::array_layout,
        typename YVector_Internal::device_type,
        typename YVector_Internal::memory_traits,
        std::is_integral<typename AMatrix_Internal::const_value_type>::value,
        false>::spmv_mv_bsrmatrix(controls, mode, alpha, A_i, x_i, beta, y_i);
    Kokkos::Profiling::popRegion();
  } else {
    Experimental::Impl::SPMV_MV_BSRMATRIX<
        typename AMatrix_Internal::const_value_type,
        typename AMatrix_Internal::const_ordinal_type,
        typename AMatrix_Internal::device_type,
        typename AMatrix_Internal::memory_traits,
        typename AMatrix_Internal::const_size_type,
        typename XVector_Internal::const_value_type**,
        typename XVector_Internal::array_layout,
        typename XVector_Internal::device_type,
        typename XVector_Internal::memory_traits,
        typename YVector_Internal::value_type**,
        typename YVector_Internal::array_layout,
        typename YVector_Internal::device_type,
        typename YVector_Internal::memory_traits,
        std::is_integral<typename AMatrix_Internal::const_value_type>::value>::
        spmv_mv_bsrmatrix(controls, mode, alpha, A_i, x_i, beta, y_i);
  }
}

/// \brief Public interface to local sparse matrix-vector multiply.
///
/// Compute y = beta*y + alpha*Op(A)*x, where x and y are either both
/// rank 1 (single vectors) or rank 2 (multivectors) Kokkos::View
/// instances, and Op(A) is determined
/// by \c mode.  If beta == 0, ignore and overwrite the initial
/// entries of y; if alpha == 0, ignore the entries of A and x.
///
/// If \c AMatrix is a KokkosSparse::Experimental::BsrMatrix, controls may have
/// \c "algorithm" = \c "experimental_bsr_tc" to use Nvidia tensor cores on
/// Volta or Ampere architectures. On Volta-architecture GPUs the only available
/// precision is mixed-precision fp32 accumulator from fp16 inputs. On
/// Ampere-architecture GPUs (cc >= 80), mixed precision is used when A is fp16,
/// x is fp16, and y is fp32. Otherwise, double-precision is used. The caller
/// may override this by setting the \c "tc_precision" = \c "mixed" or
/// \c "double" as desired.
///
/// For mixed precision, performance will degrade for blockDim < 16.
/// For double precision, for blockDim < 8.
/// For such cases, consider an alternate SpMV algorithm.
///
/// May have \c "algorithm" set to \c "native" to bypass TPLs if they are
/// enabled for Kokkos::CrsMatrix and Kokkos::Experimental::BsrMatrix on a
/// single vector, or for Kokkos::Experimental::BsrMatrix with a multivector.
///
/// \tparam AMatrix KokkosSparse::CrsMatrix or
/// KokkosSparse::Experimental::BsrMatrix
///
/// \param controls [in] kokkos-kernels control structure
/// \param mode [in] "N" for no transpose, "T" for transpose, or "C"
///   for conjugate transpose.
/// \param alpha [in] Scalar multiplier for the matrix A.
/// \param A [in] The sparse matrix A.
/// \param x [in] Either a single vector (rank-1 Kokkos::View) or
///   multivector (rank-2 Kokkos::View).
/// \param beta [in] Scalar multiplier for the (multi)vector y.
/// \param y [in/out] Either a single vector (rank-1 Kokkos::View) or
///   multivector (rank-2 Kokkos::View).  It must have the same number
///   of columns as x.
///
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
void spmv(KokkosKernels::Experimental::Controls controls, const char mode[],
          const AlphaType& alpha, const AMatrix& A, const XVector& x,
          const BetaType& beta, const YVector& y) {
  // Make sure that both x and y have the same rank.
  static_assert(
      static_cast<int>(XVector::rank) == static_cast<int>(YVector::rank),
      "KokkosSparse::spmv: Vector ranks do not match.");
  // Make sure that y is non-const.
  static_assert(std::is_same<typename YVector::value_type,
                             typename YVector::non_const_value_type>::value,
                "KokkosSparse::spmv: Output Vector must be non-const.");

  // Check compatibility of dimensions at run time.
  if ((mode[0] == NoTranspose[0]) || (mode[0] == Conjugate[0])) {
    if ((x.extent(1) != y.extent(1)) ||
        (static_cast<size_t>(A.numPointCols()) !=
         static_cast<size_t>(x.extent(0))) ||
        (static_cast<size_t>(A.numPointRows()) !=
         static_cast<size_t>(y.extent(0)))) {
      std::ostringstream os;
      os << "KokkosSparse::spmv (Generic): Dimensions do not match: "
         << ", A: " << A.numPointRows() << " x " << A.numPointCols()
         << ", x: " << x.extent(0) << " x " << x.extent(1)
         << ", y: " << y.extent(0) << " x " << y.extent(1);

      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  } else {
    if ((x.extent(1) != y.extent(1)) ||
        (static_cast<size_t>(A.numPointCols()) !=
         static_cast<size_t>(y.extent(0))) ||
        (static_cast<size_t>(A.numPointRows()) !=
         static_cast<size_t>(x.extent(0)))) {
      std::ostringstream os;
      os << "KokkosSparse::spmv (Generic): Dimensions do not match "
            "(transpose): "
         << ", A: " << A.numPointRows() << " x " << A.numPointCols()
         << ", x: " << x.extent(0) << " x " << x.extent(1)
         << ", y: " << y.extent(0) << " x " << y.extent(1);

      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  }

  if (alpha == Kokkos::ArithTraits<AlphaType>::zero() || A.numRows() == 0 ||
      A.numCols() == 0 || A.nnz() == 0) {
    // This is required to maintain semantics of KokkosKernels native SpMV:
    // if y contains NaN but beta = 0, the result y should be filled with 0.
    // For example, this is useful for passing in uninitialized y and beta=0.
    if (beta == Kokkos::ArithTraits<BetaType>::zero())
      Kokkos::deep_copy(y, Kokkos::ArithTraits<BetaType>::zero());
    else
      KokkosBlas::scal(y, beta, y);
    return;
  }
  //
  using RANK_SPECIALISE =
      typename std::conditional<static_cast<int>(XVector::rank) == 2, RANK_TWO,
                                RANK_ONE>::type;
  spmv(controls, mode, alpha, A, x, beta, y, RANK_SPECIALISE());
}

/// \brief Catch-all public interface to error on invalid Kokkos::Sparse spmv
/// argument types
///
/// This is a catch-all interfaceace that throws a compile-time error if \c
/// AMatrix is not a CrsMatrix, or BsrMatrix
///
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector,
          typename std::enable_if<
              !KokkosSparse::Experimental::is_bsr_matrix<AMatrix>::value &&
              !KokkosSparse::is_crs_matrix<AMatrix>::value>::type* = nullptr>
void spmv(KokkosKernels::Experimental::Controls /*controls*/,
          const char[] /*mode*/, const AlphaType& /*alpha*/,
          const AMatrix& /*A*/, const XVector& /*x*/, const BetaType& /*beta*/,
          const YVector& /*y*/) {
  // have to arrange this so that the compiler can't tell this is false until
  // instantiation
  static_assert(KokkosSparse::is_crs_matrix<AMatrix>::value ||
                    KokkosSparse::Experimental::is_bsr_matrix<AMatrix>::value,
                "SpMV: AMatrix must be CrsMatrix or BsrMatrix");
}

// Overload for backward compatibility and also just simpler
// interface for users that are happy with the kernel default settings
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
void spmv(const char mode[], const AlphaType& alpha, const AMatrix& A,
          const XVector& x, const BetaType& beta, const YVector& y) {
  KokkosKernels::Experimental::Controls controls;
  spmv(controls, mode, alpha, A, x, beta, y);
}

namespace Experimental {

template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
void spmv_struct(const char mode[], const int stencil_type,
                 const Kokkos::View<typename AMatrix::non_const_ordinal_type*,
                                    Kokkos::HostSpace>& structure,
                 const AlphaType& alpha, const AMatrix& A, const XVector& x,
                 const BetaType& beta, const YVector& y, const RANK_ONE) {
  // Make sure that both x and y have the same rank.
  static_assert((int)XVector::rank == (int)YVector::rank,
                "KokkosSparse::spmv_struct: Vector ranks do not match.");
  // Make sure that x (and therefore y) is rank 1.
  static_assert(
      (int)XVector::rank == 1,
      "KokkosSparse::spmv_struct: Both Vector inputs must have rank 1 in "
      "order to call this specialization of spmv.");
  // Make sure that y is non-const.
  static_assert(std::is_same<typename YVector::value_type,
                             typename YVector::non_const_value_type>::value,
                "KokkosSparse::spmv_struct: Output Vector must be non-const.");

  // Check compatibility of dimensions at run time.
  if ((mode[0] == NoTranspose[0]) || (mode[0] == Conjugate[0])) {
    if ((x.extent(1) != y.extent(1)) ||
        (static_cast<size_t>(A.numCols()) > static_cast<size_t>(x.extent(0))) ||
        (static_cast<size_t>(A.numRows()) > static_cast<size_t>(y.extent(0)))) {
      std::ostringstream os;
      os << "KokkosSparse::spmv_struct: Dimensions do not match: "
         << ", A: " << A.numRows() << " x " << A.numCols()
         << ", x: " << x.extent(0) << " x " << x.extent(1)
         << ", y: " << y.extent(0) << " x " << y.extent(1);

      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  } else {
    if ((x.extent(1) != y.extent(1)) ||
        (static_cast<size_t>(A.numCols()) > static_cast<size_t>(y.extent(0))) ||
        (static_cast<size_t>(A.numRows()) > static_cast<size_t>(x.extent(0)))) {
      std::ostringstream os;
      os << "KokkosSparse::spmv_struct: Dimensions do not match (transpose): "
         << ", A: " << A.numRows() << " x " << A.numCols()
         << ", x: " << x.extent(0) << " x " << x.extent(1)
         << ", y: " << y.extent(0) << " x " << y.extent(1);

      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  }

  typedef KokkosSparse::CrsMatrix<
      typename AMatrix::const_value_type, typename AMatrix::const_ordinal_type,
      typename AMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>,
      typename AMatrix::const_size_type>
      AMatrix_Internal;

  typedef Kokkos::View<
      typename XVector::const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<XVector>::array_layout,
      typename XVector::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
      XVector_Internal;

  typedef Kokkos::View<
      typename YVector::non_const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<YVector>::array_layout,
      typename YVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      YVector_Internal;

  AMatrix_Internal A_i = A;
  XVector_Internal x_i = x;
  YVector_Internal y_i = y;

  return KokkosSparse::Impl::SPMV_STRUCT<
      typename AMatrix_Internal::value_type,
      typename AMatrix_Internal::ordinal_type,
      typename AMatrix_Internal::device_type,
      typename AMatrix_Internal::memory_traits,
      typename AMatrix_Internal::size_type,
      typename XVector_Internal::value_type*,
      typename XVector_Internal::array_layout,
      typename XVector_Internal::device_type,
      typename XVector_Internal::memory_traits,
      typename YVector_Internal::value_type*,
      typename YVector_Internal::array_layout,
      typename YVector_Internal::device_type,
      typename YVector_Internal::memory_traits>::spmv_struct(mode, stencil_type,
                                                             structure, alpha,
                                                             A_i, x_i, beta,
                                                             y_i);
}

template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector, class XLayout = typename XVector::array_layout>
struct SPMV2D1D_STRUCT {
  static bool spmv2d1d_struct(
      const char mode[], const int stencil_type,
      const Kokkos::View<typename AMatrix::non_const_ordinal_type*,
                         Kokkos::HostSpace>& structure,
      const AlphaType& alpha, const AMatrix& A, const XVector& x,
      const BetaType& beta, const YVector& y);
};

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || !defined(KOKKOSKERNELS_ETI_ONLY)
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
struct SPMV2D1D_STRUCT<AlphaType, AMatrix, XVector, BetaType, YVector,
                       Kokkos::LayoutStride> {
  static bool spmv2d1d_struct(
      const char mode[], const int stencil_type,
      const Kokkos::View<typename AMatrix::non_const_ordinal_type*,
                         Kokkos::HostSpace>& structure,
      const AlphaType& alpha, const AMatrix& A, const XVector& x,
      const BetaType& beta, const YVector& y) {
    spmv_struct(mode, stencil_type, structure, alpha, A, x, beta, y,
                RANK_ONE());
    return true;
  }
};
#else
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
struct SPMV2D1D_STRUCT<AlphaType, AMatrix, XVector, BetaType, YVector,
                       Kokkos::LayoutStride> {
  static bool spmv2d1d_struct(
      const char /*mode*/[], const int /*stencil_type*/,
      const Kokkos::View<typename AMatrix::non_const_ordinal_type*,
                         Kokkos::HostSpace>& /*structure*/,
      const AlphaType& /*alpha*/, const AMatrix& /*A*/, const XVector& /*x*/,
      const BetaType& /*beta*/, const YVector& /*y*/) {
    return false;
  }
};
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || !defined(KOKKOSKERNELS_ETI_ONLY)
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
struct SPMV2D1D_STRUCT<AlphaType, AMatrix, XVector, BetaType, YVector,
                       Kokkos::LayoutLeft> {
  static bool spmv2d1d_struct(
      const char mode[], const int stencil_type,
      const Kokkos::View<typename AMatrix::non_const_ordinal_type*,
                         Kokkos::HostSpace>& structure,
      const AlphaType& alpha, const AMatrix& A, const XVector& x,
      const BetaType& beta, const YVector& y) {
    spmv_struct(mode, stencil_type, structure, alpha, A, x, beta, y,
                RANK_ONE());
    return true;
  }
};
#else
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
struct SPMV2D1D_STRUCT<AlphaType, AMatrix, XVector, BetaType, YVector,
                       Kokkos::LayoutLeft> {
  static bool spmv2d1d_struct(
      const char /*mode*/[], const int /*stencil_type*/,
      const Kokkos::View<typename AMatrix::non_const_ordinal_type*,
                         Kokkos::HostSpace>& /*structure*/,
      const AlphaType& /*alpha*/, const AMatrix& /*A*/, const XVector& /*x*/,
      const BetaType& /*beta*/, const YVector& /*y*/) {
    return false;
  }
};
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || !defined(KOKKOSKERNELS_ETI_ONLY)
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
struct SPMV2D1D_STRUCT<AlphaType, AMatrix, XVector, BetaType, YVector,
                       Kokkos::LayoutRight> {
  static bool spmv2d1d_struct(
      const char mode[], const int stencil_type,
      const Kokkos::View<typename AMatrix::non_const_ordinal_type*,
                         Kokkos::HostSpace>& structure,
      const AlphaType& alpha, const AMatrix& A, const XVector& x,
      const BetaType& beta, const YVector& y) {
    spmv_struct(mode, stencil_type, structure, alpha, A, x, beta, y,
                RANK_ONE());
    return true;
  }
};
#else
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
struct SPMV2D1D_STRUCT<AlphaType, AMatrix, XVector, BetaType, YVector,
                       Kokkos::LayoutRight> {
  static bool spmv2d1d_struct(
      const char /*mode*/[], const int /*stencil_type*/,
      const Kokkos::View<typename AMatrix::non_const_ordinal_type*,
                         Kokkos::HostSpace>& /*structure*/,
      const AlphaType& /*alpha*/, const AMatrix& /*A*/, const XVector& /*x*/,
      const BetaType& /*beta*/, const YVector& /*y*/) {
    return false;
  }
};
#endif

template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
void spmv_struct(const char mode[], const int stencil_type,
                 const Kokkos::View<typename AMatrix::non_const_ordinal_type*,
                                    Kokkos::HostSpace>& structure,
                 const AlphaType& alpha, const AMatrix& A, const XVector& x,
                 const BetaType& beta, const YVector& y, const RANK_TWO) {
  // Make sure that both x and y have the same rank.
  static_assert(XVector::rank == YVector::rank,
                "KokkosBlas::spmv: Vector ranks do not match.");
  // Make sure that y is non-const.
  static_assert(std::is_same<typename YVector::value_type,
                             typename YVector::non_const_value_type>::value,
                "KokkosBlas::spmv: Output Vector must be non-const.");

  // Check compatibility of dimensions at run time.
  if ((mode[0] == NoTranspose[0]) || (mode[0] == Conjugate[0])) {
    if ((x.extent(1) != y.extent(1)) ||
        (static_cast<size_t>(A.numCols()) > static_cast<size_t>(x.extent(0))) ||
        (static_cast<size_t>(A.numRows()) > static_cast<size_t>(y.extent(0)))) {
      std::ostringstream os;
      os << "KokkosBlas::spmv: Dimensions do not match: "
         << ", A: " << A.numRows() << " x " << A.numCols()
         << ", x: " << x.extent(0) << " x " << x.extent(1)
         << ", y: " << y.extent(0) << " x " << y.extent(1);
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  } else {
    if ((x.extent(1) != y.extent(1)) ||
        (static_cast<size_t>(A.numCols()) > static_cast<size_t>(y.extent(0))) ||
        (static_cast<size_t>(A.numRows()) > static_cast<size_t>(x.extent(0)))) {
      std::ostringstream os;
      os << "KokkosBlas::spmv: Dimensions do not match (transpose): "
         << ", A: " << A.numRows() << " x " << A.numCols()
         << ", x: " << x.extent(0) << " x " << x.extent(1)
         << ", y: " << y.extent(0) << " x " << y.extent(1);
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  }

  typedef KokkosSparse::CrsMatrix<
      typename AMatrix::const_value_type, typename AMatrix::const_ordinal_type,
      typename AMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>,
      typename AMatrix::const_size_type>
      AMatrix_Internal;

  AMatrix_Internal A_i = A;

  // Call single-vector version if appropriate
  if (x.extent(1) == 1) {
    typedef Kokkos::View<
        typename XVector::const_value_type*, typename YVector::array_layout,
        typename XVector::device_type,
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
        XVector_SubInternal;
    typedef Kokkos::View<
        typename YVector::non_const_value_type*, typename YVector::array_layout,
        typename YVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
        YVector_SubInternal;

    XVector_SubInternal x_i = Kokkos::subview(x, Kokkos::ALL(), 0);
    YVector_SubInternal y_i = Kokkos::subview(y, Kokkos::ALL(), 0);

    // spmv_struct (mode, alpha, A, x_i, beta, y_i);
    if (SPMV2D1D_STRUCT<AlphaType, AMatrix_Internal, XVector_SubInternal,
                        BetaType, YVector_SubInternal,
                        typename XVector_SubInternal::array_layout>::
            spmv2d1d_struct(mode, stencil_type, structure, alpha, A, x_i, beta,
                            y_i)) {
      return;
    }
  }

  // Call true rank 2 vector implementation
  {
    typedef Kokkos::View<
        typename XVector::const_value_type**, typename XVector::array_layout,
        typename XVector::device_type,
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >
        XVector_Internal;

    typedef Kokkos::View<typename YVector::non_const_value_type**,
                         typename YVector::array_layout,
                         typename YVector::device_type,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >
        YVector_Internal;

    XVector_Internal x_i = x;
    YVector_Internal y_i = y;

    return KokkosSparse::Impl::SPMV_MV<
        typename AMatrix_Internal::value_type,
        typename AMatrix_Internal::ordinal_type,
        typename AMatrix_Internal::device_type,
        typename AMatrix_Internal::memory_traits,
        typename AMatrix_Internal::size_type,
        typename XVector_Internal::value_type**,
        typename XVector_Internal::array_layout,
        typename XVector_Internal::device_type,
        typename XVector_Internal::memory_traits,
        typename YVector_Internal::value_type**,
        typename YVector_Internal::array_layout,
        typename YVector_Internal::device_type,
        typename YVector_Internal::memory_traits>::
        spmv_mv(KokkosKernels::Experimental::Controls(), mode, alpha, A_i, x_i,
                beta, y_i);
  }
}

/// \brief Public interface to structured local sparse matrix-vector multiply.
///
/// Compute y = beta*y + alpha*Op(A)*x, where x and y are either both
/// rank 1 (single vectors) or rank 2 (multivectors) Kokkos::View
/// instances, A is a KokkosSparse::CrsMatrix, and Op(A) is determined
/// by \c mode.  If beta == 0, ignore and overwrite the initial
/// entries of y; if alpha == 0, ignore the entries of A and x.
///
/// \param mode [in] "N" for no transpose, "T" for transpose, or "C"
///   for conjugate transpose.
/// \param structure [in] this 1D view stores the # rows in each dimension
/// (i,j,k) \param alpha [in] Scalar multiplier for the matrix A. \param A [in]
/// The sparse matrix; KokkosSparse::CrsMatrix instance. \param x [in] Either a
/// single vector (rank-1 Kokkos::View) or
///   multivector (rank-2 Kokkos::View).
/// \param beta [in] Scalar multiplier for the (multi)vector y.
/// \param y [in/out] Either a single vector (rank-1 Kokkos::View) or
///   multivector (rank-2 Kokkos::View).  It must have the same number
///   of columns as x.
template <class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
void spmv_struct(const char mode[], const int stencil_type,
                 const Kokkos::View<typename AMatrix::non_const_ordinal_type*,
                                    Kokkos::HostSpace>& structure,
                 const AlphaType& alpha, const AMatrix& A, const XVector& x,
                 const BetaType& beta, const YVector& y) {
  typedef
      typename std::conditional<XVector::rank == 2, RANK_TWO, RANK_ONE>::type
          RANK_SPECIALISE;
  spmv_struct(mode, stencil_type, structure, alpha, A, x, beta, y,
              RANK_SPECIALISE());
}

}  // namespace Experimental
}  // namespace KokkosSparse

#endif
