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
#ifndef KOKKOSBLAS1_ROTMG_IMPL_HPP_
#define KOKKOSBLAS1_ROTMG_IMPL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <KokkosBlas1_rotmg_spec.hpp>

namespace KokkosBlas {
namespace Impl {

template <class DXView, class YView, class PView>
KOKKOS_INLINE_FUNCTION void rotmg_impl(DXView const& d1, DXView const& d2, DXView const& x1, YView const& y1,
                                       PView const& param) {
  using Scalar = typename DXView::non_const_value_type;

  const Scalar one  = Kokkos::ArithTraits<Scalar>::one();
  const Scalar zero = Kokkos::ArithTraits<Scalar>::zero();

  const Scalar gamma      = 4096;
  const Scalar gammasq    = 4096 * 4096;
  const Scalar gammasqinv = one / gammasq;

  Scalar flag = zero;
  Scalar h11 = zero, h12 = zero, h21 = zero, h22 = zero;

  // Quick exit if d1 negative
  if (d1() < 0) {
    flag = -one;

    d1() = zero;
    d2() = zero;
    x1() = zero;
  } else {
    Scalar p2 = d2() * y1();

    // Trivial case p2 == 0
    if (p2 == zero) {
      flag     = -(one + one);
      param(0) = flag;
      return;
    }

    // General case
    Scalar p1 = d1() * x1();
    Scalar q1 = p1 * x1();
    Scalar q2 = p2 * y1();
    if (Kokkos::abs(q1) > Kokkos::abs(q2)) {
      h21      = -y1() / x1();
      h12      = p2 / p1;
      Scalar u = one - h12 * h21;
      if (u > zero) {
        flag = zero;
        d1() = d1() / u;
        d2() = d2() / u;
        x1() = x1() * u;
      } else {
        flag = -one;
        h11  = zero;
        h12  = zero;
        h21  = zero;
        h22  = zero;

        d1() = zero;
        d2() = zero;
        x1() = zero;
      }
    } else {
      if (q2 < 0) {
        flag = -one;
        h11  = zero;
        h12  = zero;
        h21  = zero;
        h22  = zero;

        d1() = zero;
        d2() = zero;
        x1() = zero;
      } else {
        flag       = one;
        h11        = p1 / p2;
        h22        = x1() / y1();
        Scalar u   = one + h11 * h22;
        Scalar tmp = d2() / u;
        d2()       = d1() / u;
        d1()       = tmp;
        x1()       = y1() * u;
      }
    }

    // Rescale d1, h11 and h12
    if (d1() != zero) {
      while ((d1() <= gammasqinv) || (d1() >= gammasq)) {
        if (flag == zero) {
          h11  = one;
          h22  = one;
          flag = -one;
        } else {
          h21  = -one;
          h12  = one;
          flag = -one;
        }

        if (d1() <= gammasqinv) {
          d1() = d1() * gammasq;
          x1() = x1() / gamma;
          h11  = h11 / gamma;
          h12  = h12 / gamma;
        } else {
          d1() = d1() / gammasq;
          x1() = x1() * gamma;
          h11  = h11 * gamma;
          h12  = h12 * gamma;
        }
      }
    }

    // Rescale d2, h21 and h22
    if (d2() != zero) {
      while ((Kokkos::abs(d2()) <= gammasqinv) || (Kokkos::abs(d2()) >= gammasq)) {
        if (flag == zero) {
          h11  = one;
          h22  = one;
          flag = -one;
        } else {
          h21  = -one;
          h12  = one;
          flag = -one;
        }

        if (Kokkos::abs(d2()) <= gammasqinv) {
          d2() = d2() * gammasq;
          h21  = h21 / gamma;
          h22  = h22 / gamma;
        } else {
          d2() = d2() / gammasq;
          h21  = h21 * gamma;
          h22  = h22 * gamma;
        }
      }
    }

    // Setup output parameters
    if (flag < zero) {
      param(1) = h11;
      param(2) = h21;
      param(3) = h12;
      param(4) = h22;
    } else if (flag == zero) {
      param(2) = h21;
      param(3) = h12;
    } else {
      param(1) = h11;
      param(4) = h22;
    }
    param(0) = flag;
  }
}

template <class DXView, class YView, class PView>
struct rotmg_functor {
  using Scalar = typename DXView::non_const_value_type;

  DXView d1, d2, x1;
  YView y1;
  PView param;

  rotmg_functor(DXView& d1_, DXView& d2_, DXView& x1_, const YView& y1_, PView& param_)
      : d1(d1_), d2(d2_), x1(x1_), y1(y1_), param(param_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int) const { rotmg_impl(d1, d2, x1, y1, param); }
};

template <class execution_space, class DXView, class YView, class PView>
void Rotmg_Invoke(execution_space const& space, DXView const& d1, DXView const& d2, DXView const& x1, YView const& y1,
                  PView const& param) {
  using Scalar = typename DXView::value_type;
  static_assert(!Kokkos::ArithTraits<Scalar>::is_complex, "rotmg is not defined for complex types!");

  rotmg_functor myFunc(d1, d2, x1, y1, param);
  Kokkos::RangePolicy<execution_space> rotmg_policy(space, 0, 1);
  Kokkos::parallel_for("KokkosBlas1::rotmg", rotmg_policy, myFunc);
}

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_ROTMG_IMPL_HPP_
