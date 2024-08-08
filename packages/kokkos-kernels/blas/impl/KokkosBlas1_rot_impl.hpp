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
#ifndef KOKKOSBLAS1_ROT_IMPL_HPP_
#define KOKKOSBLAS1_ROT_IMPL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <KokkosBlas1_rot_spec.hpp>

namespace KokkosBlas {
namespace Impl {

template <class VectorView, class ScalarView>
struct rot_functor {
  using scalar_type = typename VectorView::non_const_value_type;

  VectorView X, Y;
  ScalarView c, s;

  rot_functor(VectorView const& X_, VectorView const& Y_, ScalarView const& c_, ScalarView const& s_)
      : X(X_), Y(Y_), c(c_), s(s_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int entryIdx) const {
    const scalar_type temp = c() * X(entryIdx) + s() * Y(entryIdx);
    Y(entryIdx)            = c() * Y(entryIdx) - s() * X(entryIdx);
    X(entryIdx)            = temp;
  }
};

template <class ExecutionSpace, class VectorView, class ScalarView>
void Rot_Invoke(ExecutionSpace const& space, VectorView const& X, VectorView const& Y, ScalarView const& c,
                ScalarView const& s) {
  Kokkos::RangePolicy<ExecutionSpace> rot_policy(space, 0, X.extent(0));
  rot_functor rot_func(X, Y, c, s);
  Kokkos::parallel_for("KokkosBlas::rot", rot_policy, rot_func);
}

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_ROT_IMPL_HPP_
