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

#ifndef __KOKKOSBATCHED_TRMM_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_TRMM_SERIAL_IMPL_HPP__

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Trmm_Serial_Internal.hpp"

namespace KokkosBatched {
//// Lower non-transpose ////
template <typename ArgDiag>
struct SerialTrmm<Side::Left, Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalLeftLower<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(0), A.extent(1), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride_0(), A.stride_1(), B.data(), B.stride_0(), B.stride_1());
  }
};
template <typename ArgDiag>
struct SerialTrmm<Side::Right, Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalRightLower<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(0), A.extent(1), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride_0(), A.stride_1(), B.data(), B.stride_0(), B.stride_1());
  }
};
//// Lower transpose /////
template <typename ArgDiag>
struct SerialTrmm<Side::Left, Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalLeftUpper<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride_1(), A.stride_0(), B.data(), B.stride_0(), B.stride_1());
  }
};
template <typename ArgDiag>
struct SerialTrmm<Side::Right, Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalRightUpper<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride_1(), A.stride_0(), B.data(), B.stride_0(), B.stride_1());
  }
};
//// Lower conjugate-transpose ////
template <typename ArgDiag>
struct SerialTrmm<Side::Left, Uplo::Lower, Trans::ConjTranspose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalLeftUpper<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, true, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(), A.stride_1(),
        A.stride_0(), B.data(), B.stride_0(), B.stride_1());
  }
};
template <typename ArgDiag>
struct SerialTrmm<Side::Right, Uplo::Lower, Trans::ConjTranspose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalRightUpper<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, true, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(), A.stride_1(),
        A.stride_0(), B.data(), B.stride_0(), B.stride_1());
  }
};
//// Upper non-transpose ////
template <typename ArgDiag>
struct SerialTrmm<Side::Left, Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalLeftUpper<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(0), A.extent(1), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride_0(), A.stride_1(), B.data(), B.stride_0(), B.stride_1());
  }
};
template <typename ArgDiag>
struct SerialTrmm<Side::Right, Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalRightUpper<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(0), A.extent(1), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride_0(), A.stride_1(), B.data(), B.stride_0(), B.stride_1());
  }
};
//// Upper transpose /////
template <typename ArgDiag>
struct SerialTrmm<Side::Left, Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalLeftLower<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride_1(), A.stride_0(), B.data(), B.stride_0(), B.stride_1());
  }
};
template <typename ArgDiag>
struct SerialTrmm<Side::Right, Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalRightLower<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(),
        A.stride_1(), A.stride_0(), B.data(), B.stride_0(), B.stride_1());
  }
};
//// Upper conjugate-transpose ////
template <typename ArgDiag>
struct SerialTrmm<Side::Left, Uplo::Upper, Trans::ConjTranspose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalLeftLower<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, true, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(), A.stride_1(),
        A.stride_0(), B.data(), B.stride_0(), B.stride_1());
  }
};
template <typename ArgDiag>
struct SerialTrmm<Side::Right, Uplo::Upper, Trans::ConjTranspose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalRightLower<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, true, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(), A.stride_1(),
        A.stride_0(), B.data(), B.stride_0(), B.stride_1());
  }
};
}  // namespace KokkosBatched

#endif  // __KOKKOSBATCHED_TRMM_SERIAL_IMPL_HPP__
