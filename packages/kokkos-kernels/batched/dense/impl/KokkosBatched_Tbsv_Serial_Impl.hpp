// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_TBSV_SERIAL_IMPL_HPP_
#define KOKKOSBATCHED_TBSV_SERIAL_IMPL_HPP_

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include "KokkosBlas_util.hpp"
#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Tbsv_Serial_Internal.hpp"

namespace KokkosBatched {
namespace Impl {
template <typename AViewType, typename XViewType>
KOKKOS_INLINE_FUNCTION static int checkTbsvInput([[maybe_unused]] const AViewType &A,
                                                 [[maybe_unused]] const XViewType &x, [[maybe_unused]] const int k) {
  static_assert(Kokkos::is_view_v<AViewType>, "KokkosBatched::tbsv: AViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<XViewType>, "KokkosBatched::tbsv: XViewType is not a Kokkos::View.");
  static_assert(AViewType::rank == 2, "KokkosBatched::tbsv: AViewType must have rank 2.");
  static_assert(XViewType::rank == 1, "KokkosBatched::tbsv: XViewType must have rank 1.");

#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  if (k < 0) {
    Kokkos::printf(
        "KokkosBatched::tbsv: input parameter k must not be less than 0: k = "
        "%d\n",
        k);
    return 1;
  }

  const int lda = A.extent(0), n = A.extent(1);
  if (lda < (k + 1)) {
    Kokkos::printf(
        "KokkosBatched::tbsv: leading dimension of A must be smaller than k+1: "
        "lda = %d, k = %d\n",
        lda, k);
    return 1;
  }

  const int nx = x.extent(0);
  if (nx != n) {
    Kokkos::printf(
        "KokkosBatched::tbsv: Dimensions of x and A do not match: X: %d, A: %d "
        "x %d\n"
        "x.extent(0) must be equal to A.extent(1)\n",
        nx, lda, n);
    return 1;
  }
#endif
  return 0;
}

}  // namespace Impl

//// Lower non-transpose ////
template <typename ArgDiag>
struct SerialTbsv<Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Tbsv::Unblocked> {
  template <typename AViewType, typename XViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const XViewType &x, const int k) {
    auto info = Impl::checkTbsvInput(A, x, k);
    if (info) return info;

    return Impl::SerialTbsvInternalLower<Algo::Tbsv::Unblocked>::invoke(
        ArgDiag::use_unit_diag, A.extent(1), A.data(), A.stride(0), A.stride(1), x.data(), x.stride(0), k);
  }
};

//// Lower transpose ////
template <typename ArgDiag>
struct SerialTbsv<Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Tbsv::Unblocked> {
  template <typename AViewType, typename XViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const XViewType &x, const int k) {
    auto info = Impl::checkTbsvInput(A, x, k);
    if (info) return info;

    return Impl::SerialTbsvInternalLowerTranspose<Algo::Tbsv::Unblocked>::invoke(
        KokkosBlas::Impl::OpID(), ArgDiag::use_unit_diag, A.extent(1), A.data(), A.stride(0), A.stride(1), x.data(),
        x.stride(0), k);
  }
};

//// Lower conjugate-transpose ////
template <typename ArgDiag>
struct SerialTbsv<Uplo::Lower, Trans::ConjTranspose, ArgDiag, Algo::Tbsv::Unblocked> {
  template <typename AViewType, typename XViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const XViewType &x, const int k) {
    auto info = Impl::checkTbsvInput(A, x, k);
    if (info) return info;

    return Impl::SerialTbsvInternalLowerTranspose<Algo::Tbsv::Unblocked>::invoke(
        KokkosBlas::Impl::OpConj(), ArgDiag::use_unit_diag, A.extent(1), A.data(), A.stride(0), A.stride(1), x.data(),
        x.stride(0), k);
  }
};

//// Upper non-transpose ////
template <typename ArgDiag>
struct SerialTbsv<Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Tbsv::Unblocked> {
  template <typename AViewType, typename XViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const XViewType &x, const int k) {
    auto info = Impl::checkTbsvInput(A, x, k);
    if (info) return info;

    return Impl::SerialTbsvInternalUpper<Algo::Tbsv::Unblocked>::invoke(
        ArgDiag::use_unit_diag, A.extent(1), A.data(), A.stride(0), A.stride(1), x.data(), x.stride(0), k);
  }
};

//// Upper transpose ////
template <typename ArgDiag>
struct SerialTbsv<Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Tbsv::Unblocked> {
  template <typename AViewType, typename XViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const XViewType &x, const int k) {
    auto info = Impl::checkTbsvInput(A, x, k);
    if (info) return info;

    return Impl::SerialTbsvInternalUpperTranspose<Algo::Tbsv::Unblocked>::invoke(
        KokkosBlas::Impl::OpID(), ArgDiag::use_unit_diag, A.extent(1), A.data(), A.stride(0), A.stride(1), x.data(),
        x.stride(0), k);
  }
};

//// Upper conjugate-transpose ////
template <typename ArgDiag>
struct SerialTbsv<Uplo::Upper, Trans::ConjTranspose, ArgDiag, Algo::Tbsv::Unblocked> {
  template <typename AViewType, typename XViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const XViewType &x, const int k) {
    auto info = Impl::checkTbsvInput(A, x, k);
    if (info) return info;

    return Impl::SerialTbsvInternalUpperTranspose<Algo::Tbsv::Unblocked>::invoke(
        KokkosBlas::Impl::OpConj(), ArgDiag::use_unit_diag, A.extent(1), A.data(), A.stride(0), A.stride(1), x.data(),
        x.stride(0), k);
  }
};

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_TBSV_SERIAL_IMPL_HPP_
