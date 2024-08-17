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
#include "KokkosSparse_spmv_handle.hpp"
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

// clang-format off
/// \brief Kokkos sparse matrix-vector multiply.
/// Computes y := alpha*Op(A)*x + beta*y, where Op(A) is
/// controlled by mode (see below).
///
/// \tparam ExecutionSpace A Kokkos execution space. Must be able to access
///   the memory spaces of A, x, and y. Must match Handle::ExecutionSpaceType.
/// \tparam Handle Specialization of KokkosSparse::SPMVHandle
/// \tparam AlphaType Type of coefficient alpha. Must be convertible to
///   YVector::value_type.
/// \tparam AMatrix A KokkosSparse::CrsMatrix, or
///   KokkosSparse::Experimental::BsrMatrix. Must be identical to Handle::AMatrixType.
/// \tparam XVector Type of x, must be a rank-1 or 2 Kokkos::View. Must be identical to Handle::XVectorType.
/// \tparam BetaType Type of coefficient beta. Must be
///   convertible to YVector::value_type.
/// \tparam YVector Type of y, must be a rank-1 or 2 Kokkos::View and its rank must match that of XVector. Must
///   be identical to Handle::YVectorType.
///
/// \param space [in] The execution space instance on which to run the
///   kernel.
/// \param handle [in/out] a pointer to a KokkosSparse::SPMVHandle. On the first call to spmv with
///   a given handle instance, the handle's internal data will be initialized automatically.
///   On all later calls to spmv, this internal data will be reused.
/// \param mode [in] Select A's operator mode: "N" for normal, "T" for
///   transpose, "C" for conjugate or "H" for conjugate transpose.
/// \param alpha [in] Scalar multiplier for the matrix A.
/// \param A [in] The sparse matrix A. If handle has previously been passed to spmv, A must be identical to the
///   A passed in to that first call.
/// \param x [in] A vector to multiply on the left by A.
/// \param beta [in] Scalar multiplier for the vector y.
/// \param y [in/out] Result vector.
// clang-format on
template <class ExecutionSpace, class Handle, class AlphaType, class AMatrix, class XVector, class BetaType,
          class YVector>
void spmv(const ExecutionSpace& space, Handle* handle, const char mode[], const AlphaType& alpha, const AMatrix& A,
          const XVector& x, const BetaType& beta, const YVector& y) {
  // Make sure A is a CrsMatrix or BsrMatrix.
  static_assert(is_crs_matrix_v<AMatrix> || Experimental::is_bsr_matrix_v<AMatrix>,
                "KokkosSparse::spmv: AMatrix must be a CrsMatrix or BsrMatrix");
  // Make sure that x and y are Views.
  static_assert(Kokkos::is_view<XVector>::value, "KokkosSparse::spmv: XVector must be a Kokkos::View.");
  static_assert(Kokkos::is_view<YVector>::value, "KokkosSparse::spmv: YVector must be a Kokkos::View.");
  // Make sure A, x, y are accessible to ExecutionSpace
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename AMatrix::memory_space>::accessible,
                "KokkosSparse::spmv: AMatrix must be accessible from ExecutionSpace");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename XVector::memory_space>::accessible,
                "KokkosSparse::spmv: XVector must be accessible from ExecutionSpace");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename YVector::memory_space>::accessible,
                "KokkosSparse::spmv: YVector must be accessible from ExecutionSpace");
  // Make sure that x and y have the same rank.
  static_assert(XVector::rank() == YVector::rank(), "KokkosSparse::spmv: Vector ranks do not match.");
  // Make sure that x (and therefore y) is rank 1 or 2.
  static_assert(XVector::rank() == size_t(1) || XVector::rank() == size_t(2),
                "KokkosSparse::spmv: Both Vector inputs must have rank 1 or 2");
  // Make sure that y is non-const.
  static_assert(!std::is_const_v<typename YVector::value_type>, "KokkosSparse::spmv: Output Vector must be non-const.");

  // Check that A, X, Y types match that of the Handle
  // But only check this if Handle is the user-facing type (SPMVHandle).
  // We may internally call spmv with SPMVHandleImpl, which does not include
  // the matrix and vector types.
  if constexpr (KokkosSparse::Impl::is_spmv_handle_v<Handle>) {
    static_assert(std::is_same_v<AMatrix, typename Handle::AMatrixType>,
                  "KokkosSparse::spmv: AMatrix must be identical to Handle::AMatrixType");
    static_assert(std::is_same_v<XVector, typename Handle::XVectorType>,
                  "KokkosSparse::spmv: XVector must be identical to Handle::XVectorType");
    static_assert(std::is_same_v<YVector, typename Handle::YVectorType>,
                  "KokkosSparse::spmv: YVector must be identical to Handle::YVectorType");
  }

  constexpr bool isBSR = Experimental::is_bsr_matrix_v<AMatrix>;

  // Check compatibility of dimensions at run time.
  size_t m, n;

  if constexpr (!isBSR) {
    m = A.numRows();
    n = A.numCols();
  } else {
    m = A.numRows() * A.blockDim();
    n = A.numCols() * A.blockDim();
  }

  if ((mode[0] == NoTranspose[0]) || (mode[0] == Conjugate[0])) {
    if ((x.extent(1) != y.extent(1)) || (n != x.extent(0)) || (m != y.extent(0))) {
      std::ostringstream os;
      os << "KokkosSparse::spmv: Dimensions do not match: "
         << ", A: " << m << " x " << n << ", x: " << x.extent(0) << " x " << x.extent(1) << ", y: " << y.extent(0)
         << " x " << y.extent(1);
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  } else {
    if ((x.extent(1) != y.extent(1)) || (m != x.extent(0)) || (n != y.extent(0))) {
      std::ostringstream os;
      os << "KokkosSparse::spmv: Dimensions do not match (transpose): "
         << ", A: " << A.numRows() << " x " << A.numCols() << ", x: " << x.extent(0) << " x " << x.extent(1)
         << ", y: " << y.extent(0) << " x " << y.extent(1);
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  }

  // Efficiently handle cases where alpha*Op(A) is equivalent to the zero matrix
  if (alpha == Kokkos::ArithTraits<AlphaType>::zero() || m == 0 || n == 0 || A.nnz() == 0) {
    // This is required to maintain semantics of KokkosKernels native SpMV:
    // if y contains NaN but beta = 0, the result y should be filled with 0.
    // For example, this is useful for passing in uninitialized y and beta=0.
    if (beta == Kokkos::ArithTraits<BetaType>::zero())
      Kokkos::deep_copy(space, y, Kokkos::ArithTraits<BetaType>::zero());
    else
      KokkosBlas::scal(space, y, beta, y);
    return;
  }

  // Get the "impl" parent class of Handle, if it's not already the impl
  using HandleImpl = typename Handle::ImplType;

  using ACrs_Internal =
      CrsMatrix<typename AMatrix::const_value_type, typename AMatrix::const_ordinal_type, typename AMatrix::device_type,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>, typename AMatrix::const_size_type>;
  using ABsr_Internal =
      Experimental::BsrMatrix<typename AMatrix::const_value_type, typename AMatrix::const_ordinal_type,
                              typename AMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                              typename AMatrix::const_size_type>;

  using AMatrix_Internal = std::conditional_t<isBSR, ABsr_Internal, ACrs_Internal>;

  // Intercept special case: A is a BsrMatrix with blockDim() == 1
  // This is exactly equivalent to CrsMatrix (more performant)
  // and cuSPARSE actually errors out in that case.
  //
  // This relies on the fact that this codepath will always be taken for
  // this particular matrix (so internally, this handle is only ever used for
  // Crs)
  if constexpr (isBSR) {
    if (A.blockDim() == 1) {
      // Construct an ACrs_Internal (unmanaged memory) from A's views
      typename ACrs_Internal::row_map_type rowmap(A.graph.row_map);
      typename ACrs_Internal::index_type entries(A.graph.entries);
      typename ACrs_Internal::values_type values(A.values);
      ACrs_Internal ACrs(std::string{}, A.numRows(), A.numCols(), A.nnz(), values, rowmap, entries);
      spmv(space, handle->get_impl(), mode, alpha, ACrs, x, beta, y);
      return;
    }
  }

  AMatrix_Internal A_i(A);

  // Note: data_type of a View includes both the scalar and rank
  using XVector_Internal =
      Kokkos::View<typename XVector::const_data_type,
                   typename KokkosKernels::Impl::GetUnifiedLayout<XVector>::array_layout, typename XVector::device_type,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;

  using YVector_Internal = Kokkos::View<typename YVector::non_const_data_type,
                                        typename KokkosKernels::Impl::GetUnifiedLayout<YVector>::array_layout,
                                        typename YVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  // Special case: XVector/YVector are rank-2 but x,y both have one column and
  // are contiguous. In this case take rank-1 subviews of x,y and call the
  // rank-1 version.
  if constexpr (XVector::rank() == 2) {
    using XVector_SubInternal =
        Kokkos::View<typename XVector::const_value_type*,
                     typename KokkosKernels::Impl::GetUnifiedLayout<XVector>::array_layout,
                     typename XVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;
    using YVector_SubInternal = Kokkos::View<typename YVector::non_const_value_type*,
                                             typename KokkosKernels::Impl::GetUnifiedLayout<YVector>::array_layout,
                                             typename YVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    if (x.extent(1) == size_t(1) && x.span_is_contiguous() && y.span_is_contiguous()) {
      XVector_SubInternal xsub(x.data(), x.extent(0));
      YVector_SubInternal ysub(y.data(), y.extent(0));
      spmv(space, handle->get_impl(), mode, alpha, A, xsub, beta, ysub);
      return;
    }
  }

  XVector_Internal x_i(x);
  YVector_Internal y_i(y);

  bool useNative = is_spmv_algorithm_native(handle->get_algorithm());

  // Now call the proper implementation depending on isBSR and the rank of X/Y
  if constexpr (!isBSR) {
    if constexpr (XVector::rank() == 1) {
/////////////////
// CRS, rank 1 //
/////////////////
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
      // cuSPARSE does not support the conjugate mode (C)
      if constexpr (std::is_same_v<ExecutionSpace, Kokkos::Cuda>) {
        useNative = useNative || (mode[0] == Conjugate[0]);
      }
      // cuSPARSE 12 requires that the output (y) vector is 16-byte aligned for
      // all scalar types
#if defined(CUSPARSE_VER_MAJOR) && (CUSPARSE_VER_MAJOR == 12)
      uintptr_t yptr = uintptr_t((void*)y.data());
      if (yptr % 16 != 0) useNative = true;
#endif
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
      if (std::is_same_v<ExecutionSpace, Kokkos::HIP>) {
        useNative = useNative || (mode[0] != NoTranspose[0]);
      }
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
      if constexpr (std::is_same_v<typename AMatrix_Internal::memory_space, Kokkos::HostSpace>) {
        useNative = useNative || (mode[0] == Conjugate[0]);
      }
#ifdef KOKKOS_ENABLE_SYCL
      if constexpr (std::is_same_v<ExecutionSpace, Kokkos::Experimental::SYCL>) {
        useNative = useNative || (mode[0] != NoTranspose[0]);
      }
#endif
#endif
      if (useNative) {
        // Explicitly call the non-TPL SPMV implementation
        std::string label = "KokkosSparse::spmv[NATIVE," +
                            Kokkos::ArithTraits<typename AMatrix_Internal::non_const_value_type>::name() + "]";
        Kokkos::Profiling::pushRegion(label);
        Impl::SPMV<ExecutionSpace, HandleImpl, AMatrix_Internal, XVector_Internal, YVector_Internal, false>::spmv(
            space, handle, mode, alpha, A_i, x_i, beta, y_i);
        Kokkos::Profiling::popRegion();
      } else {
        // note: the cuSPARSE spmv wrapper defines a profiling region, so one is
        // not needed here.
        Impl::SPMV<ExecutionSpace, HandleImpl, AMatrix_Internal, XVector_Internal, YVector_Internal>::spmv(
            space, handle, mode, alpha, A_i, x_i, beta, y_i);
      }
    } else {
/////////////////
// CRS, rank 2 //
/////////////////
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
      if constexpr (std::is_same_v<ExecutionSpace, Kokkos::Cuda>) {
        useNative = useNative || (Conjugate[0] == mode[0]);
      }
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
      if constexpr (std::is_same_v<ExecutionSpace, Kokkos::HIP>) {
        useNative = useNative || (Conjugate[0] == mode[0]);
      }
#endif

      if (useNative) {
        std::string label = "KokkosSparse::spmv[NATIVE,MV," +
                            Kokkos::ArithTraits<typename AMatrix_Internal::non_const_value_type>::name() + "]";
        Kokkos::Profiling::pushRegion(label);
        return Impl::SPMV_MV<ExecutionSpace, HandleImpl, AMatrix_Internal, XVector_Internal, YVector_Internal,
                             std::is_integral<typename AMatrix_Internal::value_type>::value, false>::spmv_mv(space,
                                                                                                             handle,
                                                                                                             mode,
                                                                                                             alpha, A_i,
                                                                                                             x_i, beta,
                                                                                                             y_i);
        Kokkos::Profiling::popRegion();
      } else {
        return Impl::SPMV_MV<ExecutionSpace, HandleImpl, AMatrix_Internal, XVector_Internal, YVector_Internal>::spmv_mv(
            space, handle, mode, alpha, A_i, x_i, beta, y_i);
      }
    }
  } else {
    if constexpr (XVector::rank() == 1) {
/////////////////
// BSR, rank 1 //
/////////////////
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
      // cuSPARSE does not support the modes (C), (T), (H)
      if constexpr (std::is_same_v<ExecutionSpace, Kokkos::Cuda>) {
        useNative = useNative || (mode[0] != NoTranspose[0]);
      }
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
      if constexpr (std::is_same_v<typename AMatrix_Internal::memory_space, Kokkos::HostSpace>) {
        useNative = useNative || (mode[0] == Conjugate[0]);
      }
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
      // rocSparse does not support the modes (C), (T), (H)
      if constexpr (std::is_same_v<ExecutionSpace, Kokkos::HIP>) {
        useNative = useNative || (mode[0] != NoTranspose[0]);
      }
#endif
      if (useNative) {
        // Explicitly call the non-TPL SPMV_BSRMATRIX implementation
        std::string label = "KokkosSparse::spmv[NATIVE,BSRMATRIX," +
                            Kokkos::ArithTraits<typename AMatrix_Internal::non_const_value_type>::name() + "]";
        Kokkos::Profiling::pushRegion(label);
        Impl::SPMV_BSRMATRIX<ExecutionSpace, HandleImpl, AMatrix_Internal, XVector_Internal, YVector_Internal,
                             false>::spmv_bsrmatrix(space, handle, mode, alpha, A_i, x_i, beta, y_i);
        Kokkos::Profiling::popRegion();
      } else {
        Impl::SPMV_BSRMATRIX<ExecutionSpace, HandleImpl, AMatrix_Internal, XVector_Internal,
                             YVector_Internal>::spmv_bsrmatrix(space, handle, mode, alpha, A_i, x_i, beta, y_i);
      }
    } else {
      /////////////////
      // BSR, rank 2 //
      /////////////////
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
      // cuSPARSE does not support the modes (C), (T), (H)
      if constexpr (std::is_same_v<ExecutionSpace, Kokkos::Cuda>) {
        useNative = useNative || (mode[0] != NoTranspose[0]);
      }
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
      if constexpr (std::is_same_v<typename AMatrix_Internal::memory_space, Kokkos::HostSpace>) {
        useNative = useNative || (mode[0] == Conjugate[0]);
      }
#endif
      if (useNative) {
        // Explicitly call the non-TPL SPMV_BSRMATRIX implementation
        std::string label = "KokkosSparse::spmv[NATIVE,MV,BSMATRIX," +
                            Kokkos::ArithTraits<typename AMatrix_Internal::non_const_value_type>::name() + "]";
        Kokkos::Profiling::pushRegion(label);
        Impl::SPMV_MV_BSRMATRIX<ExecutionSpace, HandleImpl, AMatrix_Internal, XVector_Internal, YVector_Internal,
                                std::is_integral<typename AMatrix_Internal::const_value_type>::value,
                                false>::spmv_mv_bsrmatrix(space, handle, mode, alpha, A_i, x_i, beta, y_i);
        Kokkos::Profiling::popRegion();
      } else {
        Impl::SPMV_MV_BSRMATRIX<
            ExecutionSpace, HandleImpl, AMatrix_Internal, XVector_Internal, YVector_Internal,
            std::is_integral<typename AMatrix_Internal::const_value_type>::value>::spmv_mv_bsrmatrix(space, handle,
                                                                                                     mode, alpha, A_i,
                                                                                                     x_i, beta, y_i);
      }
    }
  }
}

// clang-format off
/// \brief Kokkos sparse matrix-vector multiply.
///   Computes y := alpha*Op(A)*x + beta*y, where Op(A) is controlled by mode
///   (see below).
///
/// \tparam ExecutionSpace A Kokkos execution space. Must be able to access
///   the memory spaces of A, x, and y.
/// \tparam AlphaType Type of coefficient alpha. Must be convertible to
///   YVector::value_type.
/// \tparam AMatrix A KokkosSparse::CrsMatrix, or KokkosSparse::Experimental::BsrMatrix
/// \tparam XVector Type of x, must be a rank-1 or rank-2 Kokkos::View
/// \tparam BetaType Type of coefficient beta. Must be convertible to YVector::value_type.
/// \tparam YVector Type of y, must be a Kokkos::View and its rank must match that of XVector
///
/// \param space [in] The execution space instance on which to run the kernel.
/// \param mode [in] Select A's operator mode: "N" for normal, "T" for
///   transpose, "C" for conjugate or "H" for conjugate transpose.
/// \param alpha [in] Scalar multiplier for the matrix A.
/// \param A [in] The sparse matrix A.
/// \param x [in] A vector to multiply on the left by A.
/// \param beta [in] Scalar multiplier for the vector y.
/// \param y [in/out] Result vector.
// clang-format on
template <class ExecutionSpace, class AlphaType, class AMatrix, class XVector, class BetaType, class YVector,
          typename = std::enable_if_t<Kokkos::is_execution_space<ExecutionSpace>::value>>
void spmv(const ExecutionSpace& space, const char mode[], const AlphaType& alpha, const AMatrix& A, const XVector& x,
          const BetaType& beta, const YVector& y) {
  SPMVAlgorithm algo = SPMV_FAST_SETUP;
  // Without handle reuse, native is overall faster than rocSPARSE
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
  if constexpr (std::is_same_v<typename AMatrix::execution_space, Kokkos::HIP>) algo = SPMV_NATIVE;
#endif
  SPMVHandle<typename AMatrix::execution_space, AMatrix, XVector, YVector> handle(algo);
  spmv(space, &handle, mode, alpha, A, x, beta, y);
}

// clang-format off
/// \brief Kokkos sparse matrix-vector multiply.
///   Computes y := alpha*Op(A)*x + beta*y, where Op(A) is controlled by mode
///   (see below).
///
/// \tparam Handle Specialization of KokkosSparse::SPMVHandle
/// \tparam AlphaType Type of coefficient alpha. Must be convertible to
///   YVector::value_type.
/// \tparam AMatrix A KokkosSparse::CrsMatrix, or
///   KokkosSparse::Experimental::BsrMatrix. Must be identical to Handle::AMatrixType.
/// \tparam XVector Type of x. Must be a rank-1 or 2 Kokkos::View and be identical to Handle::XVectorType.
/// \tparam BetaType Type of coefficient beta. Must be convertible to YVector::value_type.
/// \tparam YVector Type of y. Must have the same rank as XVector and be identical to Handle::YVectorType.
///
/// \param handle [in/out] a pointer to a KokkosSparse::SPMVHandle. On the first call to spmv with
///   a given handle instance, the handle's internal data will be initialized automatically.
///   On all later calls to spmv, this internal data will be reused.
/// \param mode [in] Select A's operator mode: "N" for normal, "T" for
///   transpose, "C" for conjugate or "H" for conjugate transpose.
/// \param alpha [in] Scalar multiplier for the matrix A.
/// \param A [in] The sparse matrix A.
/// \param x [in] A vector to multiply on the left by A.
/// \param beta [in] Scalar multiplier for the vector y.
/// \param y [in/out] Result vector.
// clang-format on
template <class Handle, class AlphaType, class AMatrix, class XVector, class BetaType, class YVector,
          typename = std::enable_if_t<!Kokkos::is_execution_space<Handle>::value>>
void spmv(Handle* handle, const char mode[], const AlphaType& alpha, const AMatrix& A, const XVector& x,
          const BetaType& beta, const YVector& y) {
  spmv(typename Handle::ExecutionSpaceType(), handle, mode, alpha, A, x, beta, y);
}

// clang-format off
/// \brief Kokkos sparse matrix-vector multiply.
///   Computes y := alpha*Op(A)*x + beta*y, where Op(A) is controlled by mode
///   (see below).
///
/// \tparam AlphaType Type of coefficient alpha. Must be convertible to YVector::value_type.
/// \tparam AMatrix A KokkosSparse::CrsMatrix, or KokkosSparse::Experimental::BsrMatrix
/// \tparam XVector Type of x, must be a rank-1 or rank-2 Kokkos::View
/// \tparam BetaType Type of coefficient beta. Must be convertible to YVector::value_type.
/// \tparam YVector Type of y, must be a Kokkos::View and its rank must match that of XVector
///
/// \param mode [in] Select A's operator mode: "N" for normal, "T" for
/// transpose, "C" for conjugate or "H" for conjugate transpose.
/// \param alpha [in] Scalar multiplier for the matrix A.
/// \param A [in] The sparse matrix A.
/// \param x [in] A vector to multiply on the left by A.
/// \param beta [in] Scalar multiplier for the vector y.
/// \param y [in/out] Result vector.
// clang-format on
template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
void spmv(const char mode[], const AlphaType& alpha, const AMatrix& A, const XVector& x, const BetaType& beta,
          const YVector& y) {
  SPMVAlgorithm algo = SPMV_FAST_SETUP;
  // Without handle reuse, native is overall faster than rocSPARSE
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
  if constexpr (std::is_same_v<typename AMatrix::execution_space, Kokkos::HIP>) algo = SPMV_NATIVE;
#endif
  SPMVHandle<typename AMatrix::execution_space, AMatrix, XVector, YVector> handle(algo);
  spmv(typename AMatrix::execution_space(), &handle, mode, alpha, A, x, beta, y);
}

namespace Experimental {

template <class ExecutionSpace, class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
void spmv_struct(const ExecutionSpace& space, const char mode[], const int stencil_type,
                 const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
                 const AlphaType& alpha, const AMatrix& A, const XVector& x, const BetaType& beta, const YVector& y,
                 [[maybe_unused]] const RANK_ONE& tag) {
  // Make sure that both x and y have the same rank.
  static_assert((int)XVector::rank == (int)YVector::rank, "KokkosSparse::spmv_struct: Vector ranks do not match.");
  // Make sure A, x, y are accessible to ExecutionSpace
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename AMatrix::memory_space>::accessible,
                "KokkosSparse::spmv_struct: AMatrix must be accessible from "
                "ExecutionSpace");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename XVector::memory_space>::accessible,
                "KokkosSparse::spmv_struct: XVector must be accessible from "
                "ExecutionSpace");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename YVector::memory_space>::accessible,
                "KokkosSparse::spmv_struct: YVector must be accessible from "
                "ExecutionSpace");
  // Make sure that x (and therefore y) is rank 1.
  static_assert((int)XVector::rank == 1,
                "KokkosSparse::spmv_struct: Both Vector inputs must have rank 1 in "
                "order to call this specialization of spmv.");
  // Make sure that y is non-const.
  static_assert(std::is_same_v<typename YVector::value_type, typename YVector::non_const_value_type>,
                "KokkosSparse::spmv_struct: Output Vector must be non-const.");

  // Check compatibility of dimensions at run time.
  if ((mode[0] == NoTranspose[0]) || (mode[0] == Conjugate[0])) {
    if ((x.extent(1) != y.extent(1)) || (static_cast<size_t>(A.numCols()) > static_cast<size_t>(x.extent(0))) ||
        (static_cast<size_t>(A.numRows()) > static_cast<size_t>(y.extent(0)))) {
      std::ostringstream os;
      os << "KokkosSparse::spmv_struct: Dimensions do not match: "
         << ", A: " << A.numRows() << " x " << A.numCols() << ", x: " << x.extent(0) << " x " << x.extent(1)
         << ", y: " << y.extent(0) << " x " << y.extent(1);

      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  } else {
    if ((x.extent(1) != y.extent(1)) || (static_cast<size_t>(A.numCols()) > static_cast<size_t>(y.extent(0))) ||
        (static_cast<size_t>(A.numRows()) > static_cast<size_t>(x.extent(0)))) {
      std::ostringstream os;
      os << "KokkosSparse::spmv_struct: Dimensions do not match (transpose): "
         << ", A: " << A.numRows() << " x " << A.numCols() << ", x: " << x.extent(0) << " x " << x.extent(1)
         << ", y: " << y.extent(0) << " x " << y.extent(1);

      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  }

  typedef KokkosSparse::CrsMatrix<typename AMatrix::const_value_type, typename AMatrix::const_ordinal_type,
                                  typename AMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                                  typename AMatrix::const_size_type>
      AMatrix_Internal;

  typedef Kokkos::View<typename XVector::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<XVector>::array_layout,
                       typename XVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>
      XVector_Internal;

  typedef Kokkos::View<typename YVector::non_const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<YVector>::array_layout,
                       typename YVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      YVector_Internal;

  AMatrix_Internal A_i = A;
  XVector_Internal x_i = x;
  YVector_Internal y_i = y;

  return KokkosSparse::Impl::SPMV_STRUCT<ExecutionSpace, AMatrix_Internal, XVector_Internal,
                                         YVector_Internal>::spmv_struct(space, mode, stencil_type, structure, alpha,
                                                                        A_i, x_i, beta, y_i);
}

template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
void spmv_struct(const char mode[], const int stencil_type,
                 const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
                 const AlphaType& alpha, const AMatrix& A, const XVector& x, const BetaType& beta, const YVector& y,
                 const RANK_ONE& tag) {
  spmv_struct(typename AMatrix::execution_space{}, mode, stencil_type, structure, alpha, A, x, beta, y, tag);
}

namespace Impl {
template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector,
          class XLayout = typename XVector::array_layout>
struct SPMV2D1D_STRUCT {
  static bool spmv2d1d_struct(
      const char mode[], const int stencil_type,
      const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
      const AlphaType& alpha, const AMatrix& A, const XVector& x, const BetaType& beta, const YVector& y);

  template <typename ExecutionSpace>
  static bool spmv2d1d_struct(
      const ExecutionSpace& space, const char mode[], const int stencil_type,
      const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
      const AlphaType& alpha, const AMatrix& A, const XVector& x, const BetaType& beta, const YVector& y);
};

#if defined(KOKKOSKERNELS_INST_LAYOUTSTRIDE) || !defined(KOKKOSKERNELS_ETI_ONLY)
template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
struct SPMV2D1D_STRUCT<AlphaType, AMatrix, XVector, BetaType, YVector, Kokkos::LayoutStride> {
  static bool spmv2d1d_struct(
      const char mode[], const int stencil_type,
      const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
      const AlphaType& alpha, const AMatrix& A, const XVector& x, const BetaType& beta, const YVector& y) {
    spmv_struct(mode, stencil_type, structure, alpha, A, x, beta, y, RANK_ONE());
    return true;
  }

  template <typename ExecutionSpace>
  static bool spmv2d1d_struct(
      const ExecutionSpace& space, const char mode[], const int stencil_type,
      const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
      const AlphaType& alpha, const AMatrix& A, const XVector& x, const BetaType& beta, const YVector& y) {
    spmv_struct(space, mode, stencil_type, structure, alpha, A, x, beta, y, RANK_ONE());
    return true;
  }
};
#else
template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
struct SPMV2D1D_STRUCT<AlphaType, AMatrix, XVector, BetaType, YVector, Kokkos::LayoutStride> {
  static bool spmv2d1d_struct(
      const char /*mode*/[], const int /*stencil_type*/,
      const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& /*structure*/,
      const AlphaType& /*alpha*/, const AMatrix& /*A*/, const XVector& /*x*/, const BetaType& /*beta*/,
      const YVector& /*y*/) {
    return false;
  }

  template <typename ExecutionSpace>
  static bool spmv2d1d_struct(
      const ExecutionSpace& /* space*/, const char /*mode*/[], const int /*stencil_type*/,
      const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& /*structure*/,
      const AlphaType& /*alpha*/, const AMatrix& /*A*/, const XVector& /*x*/, const BetaType& /*beta*/,
      const YVector& /*y*/) {
    return false;
  }
};
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || !defined(KOKKOSKERNELS_ETI_ONLY)
template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
struct SPMV2D1D_STRUCT<AlphaType, AMatrix, XVector, BetaType, YVector, Kokkos::LayoutLeft> {
  static bool spmv2d1d_struct(
      const char mode[], const int stencil_type,
      const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
      const AlphaType& alpha, const AMatrix& A, const XVector& x, const BetaType& beta, const YVector& y) {
    spmv_struct(mode, stencil_type, structure, alpha, A, x, beta, y, RANK_ONE());
    return true;
  }

  template <typename ExecutionSpace>
  static bool spmv2d1d_struct(
      const ExecutionSpace& space, const char mode[], const int stencil_type,
      const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
      const AlphaType& alpha, const AMatrix& A, const XVector& x, const BetaType& beta, const YVector& y) {
    spmv_struct(space, mode, stencil_type, structure, alpha, A, x, beta, y, RANK_ONE());
    return true;
  }
};
#else
template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
struct SPMV2D1D_STRUCT<AlphaType, AMatrix, XVector, BetaType, YVector, Kokkos::LayoutLeft> {
  static bool spmv2d1d_struct(
      const char /*mode*/[], const int /*stencil_type*/,
      const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& /*structure*/,
      const AlphaType& /*alpha*/, const AMatrix& /*A*/, const XVector& /*x*/, const BetaType& /*beta*/,
      const YVector& /*y*/) {
    return false;
  }

  template <typename ExecutionSpace>
  static bool spmv2d1d_struct(
      const ExecutionSpace /*space*/, const char /*mode*/[], const int /*stencil_type*/,
      const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& /*structure*/,
      const AlphaType& /*alpha*/, const AMatrix& /*A*/, const XVector& /*x*/, const BetaType& /*beta*/,
      const YVector& /*y*/) {
    return false;
  }
};
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || !defined(KOKKOSKERNELS_ETI_ONLY)
template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
struct SPMV2D1D_STRUCT<AlphaType, AMatrix, XVector, BetaType, YVector, Kokkos::LayoutRight> {
  static bool spmv2d1d_struct(
      const char mode[], const int stencil_type,
      const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
      const AlphaType& alpha, const AMatrix& A, const XVector& x, const BetaType& beta, const YVector& y) {
    spmv_struct(mode, stencil_type, structure, alpha, A, x, beta, y, RANK_ONE());
    return true;
  }

  template <typename ExecutionSpace>
  static bool spmv2d1d_struct(
      const ExecutionSpace& space, const char mode[], const int stencil_type,
      const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
      const AlphaType& alpha, const AMatrix& A, const XVector& x, const BetaType& beta, const YVector& y) {
    spmv_struct(space, mode, stencil_type, structure, alpha, A, x, beta, y, RANK_ONE());
    return true;
  }
};
#else
template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
struct SPMV2D1D_STRUCT<AlphaType, AMatrix, XVector, BetaType, YVector, Kokkos::LayoutRight> {
  static bool spmv2d1d_struct(
      const char /*mode*/[], const int /*stencil_type*/,
      const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& /*structure*/,
      const AlphaType& /*alpha*/, const AMatrix& /*A*/, const XVector& /*x*/, const BetaType& /*beta*/,
      const YVector& /*y*/) {
    return false;
  }

  template <typename ExecutionSpace>
  static bool spmv2d1d_struct(
      const ExecutionSpace& /*space*/, const char /*mode*/[], const int /*stencil_type*/,
      const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& /*structure*/,
      const AlphaType& /*alpha*/, const AMatrix& /*A*/, const XVector& /*x*/, const BetaType& /*beta*/,
      const YVector& /*y*/) {
    return false;
  }
};
#endif
}  // namespace Impl

template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector,
          class XLayout = typename XVector::array_layout>
using SPMV2D1D_STRUCT
    [[deprecated("KokkosSparse::SPMV2D1D_STRUCT is not part of the public "
                 "interface - use KokkosSparse::spmv_struct instead")]] =
        Impl::SPMV2D1D_STRUCT<AlphaType, AMatrix, XVector, BetaType, YVector>;

template <class ExecutionSpace, class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
void spmv_struct(const ExecutionSpace& space, const char mode[], const int stencil_type,
                 const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
                 const AlphaType& alpha, const AMatrix& A, const XVector& x, const BetaType& beta, const YVector& y,
                 [[maybe_unused]] const RANK_TWO& tag) {
  // Make sure A, x, y are accessible to ExecutionSpace
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename AMatrix::memory_space>::accessible,
                "KokkosSparse::spmv_struct: AMatrix must be accessible from "
                "ExecutionSpace");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename XVector::memory_space>::accessible,
                "KokkosSparse::spmv_struct: XVector must be accessible from "
                "ExecutionSpace");
  static_assert(Kokkos::SpaceAccessibility<ExecutionSpace, typename YVector::memory_space>::accessible,
                "KokkosSparse::spmv_struct: YVector must be accessible from "
                "ExecutionSpace");
  // Make sure that both x and y have the same rank.
  static_assert(XVector::rank == YVector::rank, "KokkosSparse::spmv: Vector ranks do not match.");
  // Make sure that y is non-const.
  static_assert(std::is_same_v<typename YVector::value_type, typename YVector::non_const_value_type>,
                "KokkosSparse::spmv: Output Vector must be non-const.");

  // Check compatibility of dimensions at run time.
  if ((mode[0] == NoTranspose[0]) || (mode[0] == Conjugate[0])) {
    if ((x.extent(1) != y.extent(1)) || (static_cast<size_t>(A.numCols()) > static_cast<size_t>(x.extent(0))) ||
        (static_cast<size_t>(A.numRows()) > static_cast<size_t>(y.extent(0)))) {
      std::ostringstream os;
      os << "KokkosSparse::spmv: Dimensions do not match: "
         << ", A: " << A.numRows() << " x " << A.numCols() << ", x: " << x.extent(0) << " x " << x.extent(1)
         << ", y: " << y.extent(0) << " x " << y.extent(1);
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  } else {
    if ((x.extent(1) != y.extent(1)) || (static_cast<size_t>(A.numCols()) > static_cast<size_t>(y.extent(0))) ||
        (static_cast<size_t>(A.numRows()) > static_cast<size_t>(x.extent(0)))) {
      std::ostringstream os;
      os << "KokkosSparse::spmv: Dimensions do not match (transpose): "
         << ", A: " << A.numRows() << " x " << A.numCols() << ", x: " << x.extent(0) << " x " << x.extent(1)
         << ", y: " << y.extent(0) << " x " << y.extent(1);
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  }

  typedef KokkosSparse::CrsMatrix<typename AMatrix::const_value_type, typename AMatrix::const_ordinal_type,
                                  typename AMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                                  typename AMatrix::const_size_type>
      AMatrix_Internal;

  AMatrix_Internal A_i = A;

  // Call single-vector version if appropriate
  if (x.extent(1) == 1) {
    typedef Kokkos::View<typename XVector::const_value_type*, typename YVector::array_layout,
                         typename XVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>
        XVector_SubInternal;
    typedef Kokkos::View<typename YVector::non_const_value_type*, typename YVector::array_layout,
                         typename YVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        YVector_SubInternal;

    XVector_SubInternal x_i = Kokkos::subview(x, Kokkos::ALL(), 0);
    YVector_SubInternal y_i = Kokkos::subview(y, Kokkos::ALL(), 0);

    // spmv_struct (mode, alpha, A, x_i, beta, y_i);
    if (Impl::SPMV2D1D_STRUCT<AlphaType, AMatrix_Internal, XVector_SubInternal, BetaType, YVector_SubInternal,
                              typename XVector_SubInternal::array_layout>::spmv2d1d_struct(space, mode, stencil_type,
                                                                                           structure, alpha, A, x_i,
                                                                                           beta, y_i)) {
      return;
    }
  }

  // Call true rank 2 vector implementation
  spmv(space, mode, alpha, A, x, beta, y);
}

template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
void spmv_struct(const char mode[], const int stencil_type,
                 const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
                 const AlphaType& alpha, const AMatrix& A, const XVector& x, const BetaType& beta, const YVector& y,
                 const RANK_TWO& tag) {
  spmv_struct(typename AMatrix::execution_space{}, mode, stencil_type, structure, alpha, A, x, beta, y, tag);
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
///             for conjugate transpose.
/// \param stencil_type
/// \param structure [in] this 1D view stores the # rows in each dimension
///                  (i,j,k)
/// \param alpha [in] Scalar multiplier for the matrix A.
/// \param A [in] The sparse matrix; KokkosSparse::CrsMatrix instance.
/// \param x [in] Either a
///                single vector (rank-1 Kokkos::View) or
///                multivector (rank-2 Kokkos::View).
/// \param beta [in] Scalar multiplier for the (multi)vector y.
/// \param y [in/out] Either a single vector (rank-1 Kokkos::View) or
///   multivector (rank-2 Kokkos::View).  It must have the same number
///   of columns as x.
template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
void spmv_struct(const char mode[], const int stencil_type,
                 const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
                 const AlphaType& alpha, const AMatrix& A, const XVector& x, const BetaType& beta, const YVector& y) {
  typedef typename std::conditional<XVector::rank == 2, RANK_TWO, RANK_ONE>::type RANK_SPECIALISE;
  spmv_struct(mode, stencil_type, structure, alpha, A, x, beta, y, RANK_SPECIALISE());
}

/// \brief Public interface to structured local sparse matrix-vector multiply.
///
/// Compute y = beta*y + alpha*Op(A)*x, where x and y are either both
/// rank 1 (single vectors) or rank 2 (multivectors) Kokkos::View
/// instances, A is a KokkosSparse::CrsMatrix, and Op(A) is determined
/// by \c mode.  If beta == 0, ignore and overwrite the initial
/// entries of y; if alpha == 0, ignore the entries of A and x.
///
/// \param space [in] The execution space instance on which to run the
///   kernel.
/// \param mode [in] "N" for no transpose, "T" for transpose, or "C"
///             for conjugate transpose.
/// \param stencil_type
/// \param structure [in] this 1D view stores the # rows in each dimension
///                  (i,j,k)
/// \param alpha [in] Scalar multiplier for the matrix A.
/// \param A [in] The sparse matrix; KokkosSparse::CrsMatrix instance.
/// \param x [in] Either a
///                single vector (rank-1 Kokkos::View) or
///                multivector (rank-2 Kokkos::View).
/// \param beta [in] Scalar multiplier for the (multi)vector y.
/// \param y [in/out] Either a single vector (rank-1 Kokkos::View) or
///   multivector (rank-2 Kokkos::View).  It must have the same number
///   of columns as x.
template <class ExecutionSpace, class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
void spmv_struct(const ExecutionSpace& space, const char mode[], const int stencil_type,
                 const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
                 const AlphaType& alpha, const AMatrix& A, const XVector& x, const BetaType& beta, const YVector& y) {
  typedef typename std::conditional<XVector::rank == 2, RANK_TWO, RANK_ONE>::type RANK_SPECIALISE;
  spmv_struct(space, mode, stencil_type, structure, alpha, A, x, beta, y, RANK_SPECIALISE());
}

}  // namespace Experimental
}  // namespace KokkosSparse

// Pull in all the deprecated versions of spmv
// It's included here (and not at the top) because it uses definitions in this
// file.
#include "KokkosSparse_spmv_deprecated.hpp"

#endif
