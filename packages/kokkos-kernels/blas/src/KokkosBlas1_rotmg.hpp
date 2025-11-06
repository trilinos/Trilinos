// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_ROTMG_HPP_
#define KOKKOSBLAS1_ROTMG_HPP_

#include <Kokkos_Core.hpp>
#include <KokkosBlas1_rotmg_spec.hpp>
#include <KokkosKernels_helpers.hpp>

namespace KokkosBlas {

// clang-format off
/// \brief Compute the coefficients to apply a modified Givens rotation.
///
/// \tparam execution_space the execution space where the kernel will be
///         executed
/// \tparam DXView a rank0 view type that hold non const data
/// \tparam YView a rank0 view type that holds const data
/// \tparam PView a rank1 view of
///         static extent 5 that holds non const data
///
/// \param space [in] execution space used for parallel loops
/// \param d1 [in/out] On input, square of initial x scaling factor. On output, square of x scaling factor to be applied after rotm.
/// \param d2 [in/out] On input, square of initial y scaling factor. On output, square of y scaling factor
/// to be applied after rotm.
/// \param x1 [in/out] On input, element from first vector to rotate. On output, the rotated element before scaling.
/// \param y1 [in] Element from second vector to rotate. It is not modified by this routine.
/// \param param [out] 5-element parameter array defining the rotation, to be used by rotm.
// clang-format on
template <class execution_space, class DXView, class YView, class PView>
void rotmg(execution_space const& space, DXView const& d1, DXView const& d2, DXView const& x1, YView const& y1,
           PView const& param) {
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename DXView::memory_space>::accessible,
                "rotmg: execution_space cannot access data in DXView");

  using DXView_Internal =
      Kokkos::View<typename DXView::value_type, typename KokkosKernels::Impl::GetUnifiedLayout<DXView>::array_layout,
                   Kokkos::Device<execution_space, typename DXView::memory_space>,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  using YView_Internal =
      Kokkos::View<typename YView::value_type, typename KokkosKernels::Impl::GetUnifiedLayout<YView>::array_layout,
                   Kokkos::Device<execution_space, typename YView::memory_space>,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  using PView_Internal =
      Kokkos::View<typename PView::value_type[5], typename KokkosKernels::Impl::GetUnifiedLayout<PView>::array_layout,
                   Kokkos::Device<execution_space, typename PView::memory_space>,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  DXView_Internal d1_(d1), d2_(d2), x1_(x1);
  YView_Internal y1_(y1);
  PView_Internal param_(param);

  Kokkos::Profiling::pushRegion("KokkosBlas::rotmg");
  Impl::Rotmg<execution_space, DXView_Internal, YView_Internal, PView_Internal>::rotmg(space, d1_, d2_, x1_, y1_,
                                                                                       param_);
  Kokkos::Profiling::popRegion();
}

template <class DXView, class YView, class PView>
void rotmg(DXView const& d1, DXView const& d2, DXView const& x1, YView const& y1, PView const& param) {
  const typename PView::execution_space space = typename PView::execution_space();
  rotmg(space, d1, d2, x1, y1, param);
}

}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_ROTMG_HPP_
