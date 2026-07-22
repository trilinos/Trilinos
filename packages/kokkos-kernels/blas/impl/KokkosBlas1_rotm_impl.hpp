// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBLAS1_ROTM_IMPL_HPP_
#define KOKKOSBLAS1_ROTM_IMPL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_ArithTraits.hpp>
#include <KokkosBlas1_rotm_spec.hpp>

namespace KokkosBlas {
namespace Impl {

template <class VectorView, class ParamView>
struct rotm_functor {
  using Scalar = typename VectorView::non_const_value_type;

  // Dispatch tags
  struct minus_one_tag {};
  struct zero_tag {};
  struct one_tag {};

  VectorView X, Y;
  ParamView param;

  rotm_functor(VectorView const& X_, VectorView const& Y_, ParamView const& param_) : X(X_), Y(Y_), param(param_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const minus_one_tag&, const int idx) const {
    Scalar const tmp = X(idx);
    X(idx)           = param(1) * tmp + param(3) * Y(idx);
    Y(idx)           = param(2) * tmp + param(4) * Y(idx);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const zero_tag&, const int idx) const {
    Scalar const tmp = X(idx);
    X(idx)           = tmp + param(3) * Y(idx);
    Y(idx)           = param(2) * tmp + Y(idx);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const one_tag&, const int idx) const {
    Scalar const tmp = X(idx);
    X(idx)           = param(1) * tmp + Y(idx);
    Y(idx)           = -tmp + param(4) * Y(idx);
  }
};

template <class execution_space, class VectorView, class ParamView>
void Rotm_Invoke(execution_space const& space, VectorView const& X, VectorView const& Y, ParamView const& param) {
  using Scalar = typename VectorView::value_type;
  static_assert(!KokkosKernels::ArithTraits<Scalar>::is_complex, "rotm is not defined for complex types!");

  Scalar const zero = KokkosKernels::ArithTraits<Scalar>::zero();
  Scalar const one  = KokkosKernels::ArithTraits<Scalar>::one();
  Scalar const two  = one + one;

  rotm_functor myFunc(X, Y, param);

  typename ParamView::host_mirror_type param_h = Kokkos::create_mirror_view(param);
  Kokkos::deep_copy(param_h, param);
  Scalar const flag = param_h(0);

  if (flag == -two) {
    return;
  } else if (flag == -one) {
    Kokkos::RangePolicy<execution_space, typename rotm_functor<VectorView, ParamView>::minus_one_tag> rotm_policy(
        space, 0, X.extent(0));
    Kokkos::parallel_for("KokkosBlas1::rotm_minus_one", rotm_policy, myFunc);
  } else if (flag == zero) {
    Kokkos::RangePolicy<execution_space, typename rotm_functor<VectorView, ParamView>::zero_tag> rotm_policy(
        space, 0, X.extent(0));
    Kokkos::parallel_for("KokkosBlas1::rotm_zero", rotm_policy, myFunc);
  } else if (flag == one) {
    Kokkos::RangePolicy<execution_space, typename rotm_functor<VectorView, ParamView>::one_tag> rotm_policy(
        space, 0, X.extent(0));
    Kokkos::parallel_for("KokkosBlas1::rotm_one", rotm_policy, myFunc);
  } else {
    throw std::runtime_error("KokkosBlas::rotm: param(0) is not -2, -1, 0 or 1!");
  }
}

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_ROTM_IMPL_HPP_
