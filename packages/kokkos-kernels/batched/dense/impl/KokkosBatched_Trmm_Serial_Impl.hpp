// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_TRMM_SERIAL_IMPL_HPP
#define KOKKOSBATCHED_TRMM_SERIAL_IMPL_HPP

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Trmm_Serial_Internal.hpp"

namespace KokkosBatched {
//// Lower non-transpose ////
template <typename ArgDiag>
struct SerialTrmm<Side::Left, Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalLeftLower<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(0), A.extent(1), B.extent(0), B.extent(1), alpha, A.data(), A.stride(0),
        A.stride(1), B.data(), B.stride(0), B.stride(1));
  }
};
template <typename ArgDiag>
struct SerialTrmm<Side::Right, Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalRightLower<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(0), A.extent(1), B.extent(0), B.extent(1), alpha, A.data(), A.stride(0),
        A.stride(1), B.data(), B.stride(0), B.stride(1));
  }
};
//// Lower transpose /////
template <typename ArgDiag>
struct SerialTrmm<Side::Left, Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalLeftUpper<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(), A.stride(1),
        A.stride(0), B.data(), B.stride(0), B.stride(1));
  }
};
template <typename ArgDiag>
struct SerialTrmm<Side::Right, Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalRightUpper<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(), A.stride(1),
        A.stride(0), B.data(), B.stride(0), B.stride(1));
  }
};
//// Lower conjugate-transpose ////
template <typename ArgDiag>
struct SerialTrmm<Side::Left, Uplo::Lower, Trans::ConjTranspose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalLeftUpper<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, true, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(), A.stride(1),
        A.stride(0), B.data(), B.stride(0), B.stride(1));
  }
};
template <typename ArgDiag>
struct SerialTrmm<Side::Right, Uplo::Lower, Trans::ConjTranspose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalRightUpper<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, true, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(), A.stride(1),
        A.stride(0), B.data(), B.stride(0), B.stride(1));
  }
};
//// Upper non-transpose ////
template <typename ArgDiag>
struct SerialTrmm<Side::Left, Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalLeftUpper<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(0), A.extent(1), B.extent(0), B.extent(1), alpha, A.data(), A.stride(0),
        A.stride(1), B.data(), B.stride(0), B.stride(1));
  }
};
template <typename ArgDiag>
struct SerialTrmm<Side::Right, Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalRightUpper<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(0), A.extent(1), B.extent(0), B.extent(1), alpha, A.data(), A.stride(0),
        A.stride(1), B.data(), B.stride(0), B.stride(1));
  }
};
//// Upper transpose /////
template <typename ArgDiag>
struct SerialTrmm<Side::Left, Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalLeftLower<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(), A.stride(1),
        A.stride(0), B.data(), B.stride(0), B.stride(1));
  }
};
template <typename ArgDiag>
struct SerialTrmm<Side::Right, Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalRightLower<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, false, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(), A.stride(1),
        A.stride(0), B.data(), B.stride(0), B.stride(1));
  }
};
//// Upper conjugate-transpose ////
template <typename ArgDiag>
struct SerialTrmm<Side::Left, Uplo::Upper, Trans::ConjTranspose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalLeftLower<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, true, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(), A.stride(1),
        A.stride(0), B.data(), B.stride(0), B.stride(1));
  }
};
template <typename ArgDiag>
struct SerialTrmm<Side::Right, Uplo::Upper, Trans::ConjTranspose, ArgDiag, Algo::Trmm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B) {
    return SerialTrmmInternalRightLower<Algo::Trmm::Unblocked>::invoke(
        ArgDiag::use_unit_diag, true, A.extent(1), A.extent(0), B.extent(0), B.extent(1), alpha, A.data(), A.stride(1),
        A.stride(0), B.data(), B.stride(0), B.stride(1));
  }
};
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_TRMM_SERIAL_IMPL_HPP
