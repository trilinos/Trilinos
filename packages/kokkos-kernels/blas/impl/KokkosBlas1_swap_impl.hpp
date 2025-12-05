// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
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
