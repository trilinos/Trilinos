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

#ifndef __KOKKOSBATCHED_TRTRI_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_TRTRI_SERIAL_IMPL_HPP__

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Trtri_Serial_Internal.hpp"

namespace KokkosBatched {
template <typename ArgDiag>
struct SerialTrtri<Uplo::Lower, ArgDiag, Algo::Trtri::Unblocked> {
  template <typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A) {
    return SerialTrtriInternalLower<Algo::Trtri::Unblocked>::invoke(ArgDiag::use_unit_diag, A.extent(0), A.extent(1),
                                                                    A.data(), A.stride_0(), A.stride_1());
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

#endif  // __KOKKOSBATCHED_TRTRI_SERIAL_IMPL_HPP__
