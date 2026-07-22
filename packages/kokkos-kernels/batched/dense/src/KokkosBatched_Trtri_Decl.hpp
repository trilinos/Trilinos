// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_TRTRI_DECL_HPP
#define KOKKOSBATCHED_TRTRI_DECL_HPP

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

template <typename ArgUplo, typename ArgDiag, typename ArgAlgo>
struct SerialTrtri {
  template <typename ScalarType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A);
};
}  // namespace KokkosBatched
#endif  // KOKKOSBATCHED_TRTRI_DECL_HPP
