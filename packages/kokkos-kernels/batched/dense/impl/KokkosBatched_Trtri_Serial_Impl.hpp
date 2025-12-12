// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_TRTRI_SERIAL_IMPL_HPP
#define KOKKOSBATCHED_TRTRI_SERIAL_IMPL_HPP

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Trtri_Serial_Internal.hpp"

namespace KokkosBatched {
template <typename ArgDiag>
struct SerialTrtri<Uplo::Lower, ArgDiag, Algo::Trtri::Unblocked> {
  template <typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A) {
    return SerialTrtriInternalLower<Algo::Trtri::Unblocked>::invoke(ArgDiag::use_unit_diag, A.extent(0), A.extent(1),
                                                                    A.data(), A.stride(0), A.stride(1));
  }
};
template <typename ArgDiag>
struct SerialTrtri<Uplo::Upper, ArgDiag, Algo::Trtri::Unblocked> {
  template <typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A) {
    return SerialTrtriInternalUpper<Algo::Trtri::Unblocked>::invoke(ArgDiag::use_unit_diag, A.extent(0), A.extent(1),
                                                                    A.data(), A.stride(0), A.stride(1));
  }
};
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_TRTRI_SERIAL_IMPL_HPP
