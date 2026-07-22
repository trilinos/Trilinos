// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_ROTG_HPP_
#define KOKKOSBLAS1_ROTG_HPP_

#include <Kokkos_Core.hpp>
#include <KokkosBlas1_rotg_spec.hpp>
#include <KokkosKernels_helpers.hpp>

namespace KokkosBlas {

/// \brief Compute the coefficients to apply a Givens rotation.
///
/// \tparam execution_space Space on which to execute. Must be able to access SViewType and MViewType.
/// \tparam SViewType 0-D View type containing a nonconst real or complex scalar
/// \tparam MViewType 0-D View type containing a nonconst real/magnitude-typed scalar
///
/// \param space [in] the execution space
/// \param a [in/out] on input one of the values to rotate, on output the
///          rotated value r
/// \param b [in/out] on input one of the values to rotate, on
///          output the rotation parameter z
/// \param c [out] cosine value associated with the
///          rotation
/// \param s [out] sine value associated with the rotation
template <class execution_space, class SViewType, class MViewType>
void rotg(execution_space const& space, SViewType const& a, SViewType const& b, MViewType const& c,
          SViewType const& s) {
  static_assert(SViewType::rank == 0, "rotg: the inputs need to be rank 0 views");
  static_assert(MViewType::rank == 0, "rotg: the inputs need to be rank 0 views");
  static_assert(!KokkosKernels::ArithTraits<typename MViewType::value_type>::is_complex);
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename SViewType::memory_space>::accessible,
                "rotg: execution_space cannot access data in SViewType");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename MViewType::memory_space>::accessible,
                "rotg: execution_space cannot access data in MViewType");
  static_assert(!KokkosKernels::ArithTraits<typename MViewType::value_type>::is_complex,
                "rotg: MViewType cannot hold complex values.");

  using SView_Internal = Kokkos::View<
      typename SViewType::value_type, typename KokkosKernels::Impl::GetUnifiedLayout<SViewType>::array_layout,
      Kokkos::Device<execution_space, typename SViewType::memory_space>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using MView_Internal = Kokkos::View<
      typename MViewType::value_type, typename KokkosKernels::Impl::GetUnifiedLayout<MViewType>::array_layout,
      Kokkos::Device<execution_space, typename MViewType::memory_space>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  SView_Internal a_(a), b_(b), s_(s);
  MView_Internal c_(c);

  Kokkos::Profiling::pushRegion("KokkosBlas::rotg");
  Impl::Rotg<execution_space, SView_Internal, MView_Internal>::rotg(space, a, b, c, s);
  Kokkos::Profiling::popRegion();
}

}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_ROTG_HPP_
