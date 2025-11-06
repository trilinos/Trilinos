// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBLAS1_ROT_IMPL_HPP_
#define KOKKOSBLAS1_ROT_IMPL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <KokkosBlas1_rot_spec.hpp>

namespace KokkosBlas {
namespace Impl {

template <class VectorView, class MagnitudeView, class ScalarView>
struct rot_functor {
  using scalar_type = typename VectorView::non_const_value_type;

  VectorView X, Y;
  MagnitudeView c;
  ScalarView s;

  rot_functor(VectorView const& X_, VectorView const& Y_, MagnitudeView const& c_, ScalarView const& s_)
      : X(X_), Y(Y_), c(c_), s(s_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int entryIdx) const {
    const scalar_type temp = c() * X(entryIdx) + s() * Y(entryIdx);
    Y(entryIdx)            = c() * Y(entryIdx) - s() * X(entryIdx);
    X(entryIdx)            = temp;
  }
};

template <class ExecutionSpace, class VectorView, class MagnitudeView, class ScalarView>
void Rot_Invoke(ExecutionSpace const& space, VectorView const& X, VectorView const& Y, MagnitudeView const& c,
                ScalarView const& s) {
  Kokkos::RangePolicy<ExecutionSpace> rot_policy(space, 0, X.extent(0));
  rot_functor rot_func(X, Y, c, s);
  Kokkos::parallel_for("KokkosBlas::rot", rot_policy, rot_func);
}

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_ROT_IMPL_HPP_
