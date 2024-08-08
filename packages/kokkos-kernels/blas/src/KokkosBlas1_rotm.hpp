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

#ifndef KOKKOSBLAS1_ROTM_HPP_
#define KOKKOSBLAS1_ROTM_HPP_

#include <Kokkos_Core.hpp>
#include <KokkosBlas1_rotm_spec.hpp>

namespace KokkosBlas {

/// \brief Applies modified Givens rotation coefficients to vectors x and y.
///
/// \tparam execution_space the execution space where the kernel will be
///         executed, it can be used to specify a stream too.
/// \tparam VectorView a rank1 view type that hold non const data
/// \tparam ParamView a rank1 view of static extent [5] type that
///         holds const data
///
/// \param space [in]  execution space used for parallel loops in this kernel
/// \param X [in/out] vector to be rotated with param coefficients
/// \param Y [in/out] vector to be rotated with param coefficients
/// \param param [in]  output of rotmg contains rotation coefficients
///
template <class execution_space, class VectorView, class ParamView>
void rotm(execution_space const& space, VectorView const& X, VectorView const& Y, ParamView const& param) {
  static_assert(Kokkos::is_execution_space<execution_space>::value,
                "rotm: execution_space template parameter is not a Kokkos "
                "execution space.");
  static_assert(VectorView::rank == 1, "rotm: VectorView template parameter needs to be a rank 1 view");
  static_assert(ParamView::rank == 1, "rotm: ParamView template parameter needs to be a rank 1 view");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename VectorView::memory_space>::accessible,
                "rotm: VectorView template parameter memory space needs to be accessible "
                "from execution_space template parameter");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename ParamView::memory_space>::accessible,
                "rotm: ScalarView template parameter memory space needs to be accessible "
                "from execution_space template parameter");
  static_assert(std::is_same<typename VectorView::non_const_value_type, typename VectorView::value_type>::value,
                "rotm: VectorView template parameter needs to store non-const values");
  static_assert(!Kokkos::ArithTraits<typename VectorView::value_type>::is_complex,
                "rotm: VectorView template parameter cannot use complex value_type");
  static_assert(!Kokkos::ArithTraits<typename ParamView::value_type>::is_complex,
                "rotm: ParamView template parameter cannot use complex value_type");

  using VectorView_Internal = Kokkos::View<typename VectorView::non_const_value_type*,
                                           typename KokkosKernels::Impl::GetUnifiedLayout<VectorView>::array_layout,
                                           Kokkos::Device<execution_space, typename VectorView::memory_space>,
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  using ParamView_Internal = Kokkos::View<
      typename ParamView::const_value_type[5], typename KokkosKernels::Impl::GetUnifiedLayout<ParamView>::array_layout,
      Kokkos::Device<execution_space, typename ParamView::memory_space>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  VectorView_Internal X_(X), Y_(Y);
  ParamView_Internal param_(param);

  Kokkos::Profiling::pushRegion("KokkosBlas::rotm");
  Impl::Rotm<execution_space, VectorView_Internal, ParamView_Internal>::rotm(space, X_, Y_, param_);
  Kokkos::Profiling::popRegion();
}

template <class VectorView, class ParamView>
void rotm(VectorView const& X, VectorView const& Y, ParamView const& param) {
  const typename VectorView::execution_space space = typename VectorView::execution_space();
  rotm(space, X, Y, param);
}

}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_ROTM_HPP_
