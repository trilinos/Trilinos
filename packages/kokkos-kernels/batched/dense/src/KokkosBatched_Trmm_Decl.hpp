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

#ifndef KOKKOSBATCHED_TRMM_DECL_HPP
#define KOKKOSBATCHED_TRMM_DECL_HPP

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

template <typename ArgSide, typename ArgUplo, typename ArgTrans, typename ArgDiag, typename ArgAlgo>
struct SerialTrmm {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B);
};
}  // namespace KokkosBatched
#endif  // KOKKOSBATCHED_TRMM_DECL_HPP
