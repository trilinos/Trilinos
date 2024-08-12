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

#ifndef KOKKOSBATCHED_TBSV_SERIAL_IMPL_HPP_
#define KOKKOSBATCHED_TBSV_SERIAL_IMPL_HPP_

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Tbsv_Serial_Internal.hpp"

namespace KokkosBatched {

template <typename AViewType, typename XViewType>
KOKKOS_INLINE_FUNCTION static int checkTbsvInput([[maybe_unused]] const AViewType &A,
                                                 [[maybe_unused]] const XViewType &x, [[maybe_unused]] const int k) {
  static_assert(Kokkos::is_view<AViewType>::value, "KokkosBatched::tbsv: AViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view<XViewType>::value, "KokkosBatched::tbsv: XViewType is not a Kokkos::View.");
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

//// Lower non-transpose ////
template <typename ArgDiag>
struct SerialTbsv<Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Tbsv::Unblocked> {
  template <typename AViewType, typename XViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const XViewType &x, const int k) {
    auto info = checkTbsvInput(A, x, k);
    if (info) return info;

    return SerialTbsvInternalLower<Algo::Tbsv::Unblocked>::invoke(
        ArgDiag::use_unit_diag, A.extent(1), A.data(), A.stride_0(), A.stride_1(), x.data(), x.stride_0(), k);
  }
};

//// Lower transpose ////
template <typename ArgDiag>
struct SerialTbsv<Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Tbsv::Unblocked> {
  template <typename AViewType, typename XViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const XViewType &x, const int k) {
    auto info = checkTbsvInput(A, x, k);
    if (info) return info;

    return SerialTbsvInternalLowerTranspose<Algo::Tbsv::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(1), A.data(), A.stride_0(), A.stride_1(), x.data(), x.stride_0(), k);
  }
};

//// Lower conjugate-transpose ////
template <typename ArgDiag>
struct SerialTbsv<Uplo::Lower, Trans::ConjTranspose, ArgDiag, Algo::Tbsv::Unblocked> {
  template <typename AViewType, typename XViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const XViewType &x, const int k) {
    auto info = checkTbsvInput(A, x, k);
    if (info) return info;

    return SerialTbsvInternalLowerTranspose<Algo::Tbsv::Unblocked>::invoke(
        ArgDiag::use_unit_diag, true, A.extent(1), A.data(), A.stride_0(), A.stride_1(), x.data(), x.stride_0(), k);
  }
};

//// Upper non-transpose ////
template <typename ArgDiag>
struct SerialTbsv<Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Tbsv::Unblocked> {
  template <typename AViewType, typename XViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const XViewType &x, const int k) {
    auto info = checkTbsvInput(A, x, k);
    if (info) return info;

    return SerialTbsvInternalUpper<Algo::Tbsv::Unblocked>::invoke(
        ArgDiag::use_unit_diag, A.extent(1), A.data(), A.stride_0(), A.stride_1(), x.data(), x.stride_0(), k);
  }
};

//// Upper transpose ////
template <typename ArgDiag>
struct SerialTbsv<Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Tbsv::Unblocked> {
  template <typename AViewType, typename XViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const XViewType &x, const int k) {
    auto info = checkTbsvInput(A, x, k);
    if (info) return info;

    return SerialTbsvInternalUpperTranspose<Algo::Tbsv::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(1), A.data(), A.stride_0(), A.stride_1(), x.data(), x.stride_0(), k);
  }
};

//// Upper conjugate-transpose ////
template <typename ArgDiag>
struct SerialTbsv<Uplo::Upper, Trans::ConjTranspose, ArgDiag, Algo::Tbsv::Unblocked> {
  template <typename AViewType, typename XViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const XViewType &x, const int k) {
    auto info = checkTbsvInput(A, x, k);
    if (info) return info;

    return SerialTbsvInternalUpperTranspose<Algo::Tbsv::Unblocked>::invoke(
        ArgDiag::use_unit_diag, true, A.extent(1), A.data(), A.stride_0(), A.stride_1(), x.data(), x.stride_0(), k);
  }
};

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_TBSV_SERIAL_IMPL_HPP_
