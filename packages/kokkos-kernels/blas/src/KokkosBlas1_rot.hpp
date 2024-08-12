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

#ifndef KOKKOSBLAS1_ROT_HPP_
#define KOKKOSBLAS1_ROT_HPP_

#include <KokkosBlas1_rot_spec.hpp>

namespace KokkosBlas {

template <class execution_space, class VectorView, class ScalarView>
void rot(execution_space const& space, VectorView const& X, VectorView const& Y, ScalarView const& c,
         ScalarView const& s) {
  static_assert(Kokkos::is_execution_space<execution_space>::value,
                "rot: execution_space template parameter is not a Kokkos "
                "execution space.");
  static_assert(VectorView::rank == 1, "rot: VectorView template parameter needs to be a rank 1 view");
  static_assert(ScalarView::rank == 0, "rot: ScalarView template parameter needs to be a rank 0 view");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename VectorView::memory_space>::accessible,
                "rot: VectorView template parameter memory space needs to be accessible "
                "from "
                "execution_space template parameter");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename ScalarView::memory_space>::accessible,
                "rot: VectorView template parameter memory space needs to be accessible "
                "from "
                "execution_space template parameter");
  static_assert(std::is_same<typename VectorView::non_const_value_type, typename VectorView::value_type>::value,
                "rot: VectorView template parameter needs to store non-const values");

  using VectorView_Internal = Kokkos::View<typename VectorView::non_const_value_type*,
                                           typename KokkosKernels::Impl::GetUnifiedLayout<VectorView>::array_layout,
                                           Kokkos::Device<execution_space, typename VectorView::memory_space>,
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  using ScalarView_Internal = Kokkos::View<typename ScalarView::non_const_value_type,
                                           typename KokkosKernels::Impl::GetUnifiedLayout<ScalarView>::array_layout,
                                           Kokkos::Device<execution_space, typename ScalarView::memory_space>,
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  VectorView_Internal X_(X), Y_(Y);
  ScalarView_Internal c_(c), s_(s);

  Kokkos::Profiling::pushRegion("KokkosBlas::rot");
  Impl::Rot<execution_space, VectorView_Internal, ScalarView_Internal>::rot(space, X_, Y_, c_, s_);
  Kokkos::Profiling::popRegion();
}

template <class VectorView, class ScalarView>
void rot(VectorView const& X, VectorView const& Y, ScalarView const& c, ScalarView const& s) {
  const typename VectorView::execution_space space = typename VectorView::execution_space();
  rot(space, X, Y, c, s);
}

}  // namespace KokkosBlas
#endif  // KOKKOSBLAS1_ROT_HPP_
