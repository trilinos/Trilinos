/*
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
*/
#ifndef KOKKOSBLAS1_SWAP_IMPL_HPP_
#define KOKKOSBLAS1_SWAP_IMPL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>

namespace KokkosBlas {
namespace Impl {

template <class XVector, class YVector>
struct swap_functor {
  using scalar_type = typename XVector::non_const_value_type;

  XVector X;
  YVector Y;

  swap_functor(XVector const& X_, YVector const& Y_) : X(X_), Y(Y_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int const entryIdx) const {
    scalar_type const temp = Y(entryIdx);
    Y(entryIdx)            = X(entryIdx);
    X(entryIdx)            = temp;
  }
};

template <class ExecutionSpace, class XVector, class YVector>
void Swap_Invoke(ExecutionSpace const& space, XVector const& X, YVector const& Y) {
  Kokkos::RangePolicy<ExecutionSpace> swap_policy(space, 0, X.extent(0));
  swap_functor swap_func(X, Y);
  Kokkos::parallel_for("KokkosBlas::swap", swap_policy, swap_func);
}

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_SWAP_IMPL_HPP_
