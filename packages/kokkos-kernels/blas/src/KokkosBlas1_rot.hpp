// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_ROT_HPP_
#define KOKKOSBLAS1_ROT_HPP_

#include <Kokkos_Core.hpp>
#include <KokkosBlas1_rot_spec.hpp>
#include <KokkosKernels_helpers.hpp>
#include <KokkosKernels_Error.hpp>

namespace KokkosBlas {

// clang-format off
/// \brief Apply a plane rotation.
///
/// \tparam execution_space Space on which to execute. Must be able to access VectorView, MagnitudeView and ScalarView.
/// \tparam VectorView A 1-D view of nonconst scalars
/// \tparam MagnitudeView A 0-D view of nonconst, real-valued scalar
/// \tparam ScalarView A 0-D view of scalar
///
/// \param space [in] the execution space
/// \param X [in/out] First vector to be rotated
/// \param Y [in/out] Second vector to be rotated
/// \param c [out] cosine value associated with the
///          rotation
/// \param s [out] sine value associated with the rotation
// clang-format on
template <class execution_space, class VectorView, class MagnitudeView, class ScalarView>
void rot(execution_space const& space, VectorView const& X, VectorView const& Y, MagnitudeView const& c,
         ScalarView const& s) {
  static_assert(Kokkos::is_execution_space<execution_space>::value,
                "rot: execution_space template parameter is not a Kokkos "
                "execution space.");
  static_assert(Kokkos::is_view_v<VectorView>, "KokkosBlas::rot: VectorView is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<MagnitudeView>, "KokkosBlas::rot: MagnitudeView is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<ScalarView>, "KokkosBlas::rot: ScalarView is not a Kokkos::View.");
  static_assert(VectorView::rank == 1, "rot: VectorView template parameter needs to be a rank 1 view");
  static_assert(MagnitudeView::rank == 0, "rot: MagnitudeView template parameter needs to be a rank 0 view");
  static_assert(ScalarView::rank == 0, "rot: ScalarView template parameter needs to be a rank 0 view");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename VectorView::memory_space>::accessible,
                "rot: VectorView template parameter memory space needs to be accessible "
                "from "
                "execution_space template parameter");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename MagnitudeView::memory_space>::accessible,
                "rot: MagnitudeView template parameter memory space needs to be accessible "
                "from "
                "execution_space template parameter");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename ScalarView::memory_space>::accessible,
                "rot: ScalarView template parameter memory space needs to be accessible "
                "from "
                "execution_space template parameter");
  static_assert(std::is_same<typename VectorView::non_const_value_type, typename VectorView::value_type>::value,
                "rot: VectorView template parameter needs to store non-const values");

  // Check compatibility of dimensions at run time.
  if (X.extent(0) != Y.extent(0)) {
    std::ostringstream os;
    os << "KokkosBlas::rot: Dimensions of X and Y do not match: "
       << "X: " << X.extent(0) << ", Y: " << Y.extent(0);
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  using VectorView_Internal = Kokkos::View<typename VectorView::non_const_value_type*,
                                           typename KokkosKernels::Impl::GetUnifiedLayout<VectorView>::array_layout,
                                           Kokkos::Device<execution_space, typename VectorView::memory_space>,
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  using MagnitudeView_Internal = Kokkos::View<typename MagnitudeView::non_const_value_type,
                                              typename KokkosKernels::Impl::GetUnifiedLayout<ScalarView>::array_layout,
                                              Kokkos::Device<execution_space, typename ScalarView::memory_space>,
                                              Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  using ScalarView_Internal = Kokkos::View<typename ScalarView::non_const_value_type,
                                           typename KokkosKernels::Impl::GetUnifiedLayout<ScalarView>::array_layout,
                                           Kokkos::Device<execution_space, typename ScalarView::memory_space>,
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  VectorView_Internal X_(X), Y_(Y);
  MagnitudeView_Internal c_(c);
  ScalarView_Internal s_(s);

  Kokkos::Profiling::pushRegion("KokkosBlas::rot");
  Impl::Rot<execution_space, VectorView_Internal, MagnitudeView_Internal, ScalarView_Internal>::rot(space, X_, Y_, c_,
                                                                                                    s_);
  Kokkos::Profiling::popRegion();
}

template <class VectorView, class MagnitudeView, class ScalarView>
void rot(VectorView const& X, VectorView const& Y, MagnitudeView const& c, ScalarView const& s) {
  const typename VectorView::execution_space space = typename VectorView::execution_space();
  rot(space, X, Y, c, s);
}

}  // namespace KokkosBlas
#endif  // KOKKOSBLAS1_ROT_HPP_
