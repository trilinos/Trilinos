// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_GETRS_SERIAL_INTERNAL_HPP_
#define KOKKOSBATCHED_GETRS_SERIAL_INTERNAL_HPP_

#include <KokkosBatched_Util.hpp>

namespace KokkosBatched {
namespace Impl {

template <typename ArgTrans, typename ArgAlgo>
struct SerialGetrsInternal {
  template <typename AViewType, typename PivViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const PivViewType &piv, const BViewType &b);
};

//// Non-transpose ////
template <>
template <typename AViewType, typename PivViewType, typename BViewType>
KOKKOS_INLINE_FUNCTION int SerialGetrsInternal<Trans::NoTranspose, Algo::Getrs::Unblocked>::invoke(
    const AViewType &A, const PivViewType &piv, const BViewType &b) {
  KokkosBatched::SerialLaswp<Direct::Forward>::invoke(piv, b);
  KokkosBatched::SerialTrsm<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit, Algo::Trsm::Unblocked>::invoke(
      1.0, A, b);
  KokkosBatched::SerialTrsm<Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit, Algo::Trsm::Unblocked>::invoke(
      1.0, A, b);

  return 0;
}

//// Transpose ////
template <>
template <typename AViewType, typename PivViewType, typename BViewType>
KOKKOS_INLINE_FUNCTION int SerialGetrsInternal<Trans::Transpose, Algo::Getrs::Unblocked>::invoke(const AViewType &A,
                                                                                                 const PivViewType &piv,
                                                                                                 const BViewType &b) {
  KokkosBatched::SerialTrsm<Side::Left, Uplo::Upper, Trans::Transpose, Diag::NonUnit, Algo::Trsm::Unblocked>::invoke(
      1.0, A, b);
  KokkosBatched::SerialTrsm<Side::Left, Uplo::Lower, Trans::Transpose, Diag::Unit, Algo::Trsm::Unblocked>::invoke(1.0,
                                                                                                                  A, b);
  KokkosBatched::SerialLaswp<Direct::Backward>::invoke(piv, b);

  return 0;
}

//// Conj-Transpose ////
template <>
template <typename AViewType, typename PivViewType, typename BViewType>
KOKKOS_INLINE_FUNCTION int SerialGetrsInternal<Trans::ConjTranspose, Algo::Getrs::Unblocked>::invoke(
    const AViewType &A, const PivViewType &piv, const BViewType &b) {
  KokkosBatched::SerialTrsm<Side::Left, Uplo::Upper, Trans::ConjTranspose, Diag::NonUnit,
                            Algo::Trsm::Unblocked>::invoke(1.0, A, b);
  KokkosBatched::SerialTrsm<Side::Left, Uplo::Lower, Trans::ConjTranspose, Diag::Unit, Algo::Trsm::Unblocked>::invoke(
      1.0, A, b);
  KokkosBatched::SerialLaswp<Direct::Backward>::invoke(piv, b);

  return 0;
}

}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_GETRS_SERIAL_INTERNAL_HPP_
