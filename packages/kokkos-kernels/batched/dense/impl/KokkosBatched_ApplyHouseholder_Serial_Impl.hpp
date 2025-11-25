// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_APPLY_HOUSEHOLDER_SERIAL_IMPL_HPP
#define KOKKOSBATCHED_APPLY_HOUSEHOLDER_SERIAL_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)
/// \author Luc Berger-Vergiat (lberge@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Householder_Serial_Internal.hpp"

namespace KokkosBatched {

///
/// Serial Impl
/// ===========

template <typename ArgTrans>
struct SerialApplyHouseholder<Side::Left, ArgTrans> {
  template <typename uViewType, typename tauViewType, typename AViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const uViewType &u2, const tauViewType &tau, const AViewType &A,
                                           const wViewType &w) {
    if constexpr (AViewType::rank() == 1) {
      return SerialApplyLeftHouseholderInternal<ArgTrans>::invoke(A.extent(0) - 1, 1, tau.data(), u2.data(),
                                                                  u2.stride(0), A.data(), 1, A.data() + A.stride(0),
                                                                  A.stride(0), 1, w.data());
    } else {
      return SerialApplyLeftHouseholderInternal<ArgTrans>::invoke(
          A.extent(0) - 1, A.extent(1), tau.data(), u2.data(), u2.stride(0), A.data(), A.stride(1),
          A.data() + A.stride(0), A.stride(0), A.stride(1), w.data());
    }
  }
};

template <typename ArgTrans>
struct SerialApplyHouseholder<Side::Right, ArgTrans> {
  template <typename uViewType, typename tauViewType, typename AViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const uViewType &u2, const tauViewType &tau, const AViewType &A,
                                           const wViewType &w) {
    return SerialApplyRightHouseholderInternal::invoke(A.extent(0), A.extent(1) - 1, tau.data(), u2.data(),
                                                       u2.stride(0), A.data(), A.stride(0), A.data() + A.stride(1),
                                                       A.stride(0), A.stride(1), w.data());
  }
};

}  // namespace KokkosBatched

#endif
