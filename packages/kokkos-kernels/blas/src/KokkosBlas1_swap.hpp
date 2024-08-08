/*
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
*/

#ifndef KOKKOSBLAS1_SWAP_HPP_
#define KOKKOSBLAS1_SWAP_HPP_

#include <KokkosBlas1_swap_spec.hpp>

namespace KokkosBlas {

/// \brief Swaps the entries of vectors x and y.
///
/// \tparam execution_space an execution space to perform parallel work
/// \tparam XVector Type of the first vector x; a rank 1 Kokkos::View.
/// \tparam YVector Type of the first vector y; a rank 1 Kokkos::View.
///
/// \param space [in] execution space passed to execution policies
/// \param x [in/out] rank 1 View.
/// \param y [in/out] rank 1 View.
///
/// Swaps x and y. Note that this is akin to performing a deep_copy, swapping
/// pointers inside view can only be performed if no aliasing, subviews, etc...
/// exist, which cannot be asserted by this function.
///
/// This function is non-blocking unless the underlying TPL requested
/// at compile time is itself blocking
template <class execution_space, class XVector, class YVector>
void swap(execution_space const& space, XVector const& x, YVector const& y) {
  // Assert properties of XVector
  static_assert(Kokkos::is_view<XVector>::value, "KokkosBlas::swap: XVector must be a Kokkos::View.");
  static_assert(XVector::rank == 1,
                "KokkosBlas::swap: "
                "Input vector x must have rank 1.");
  static_assert(std::is_same_v<typename XVector::value_type, typename XVector::non_const_value_type>,
                "KokkosBlas::swap: XVector must store non const values.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename XVector::memory_space>::accessible,
                "swap: execution_space cannot access data in XVector");

  // Assert properties of YVector, could probably use a function for this as
  // XVector and YVector are required to have identical properties...
  static_assert(Kokkos::is_view<YVector>::value, "KokkosBlas::swap: YVector must be a Kokkos::View.");
  static_assert(YVector::rank == 1,
                "KokkosBlas::swap: "
                "Input vector y must have rank 1.");
  static_assert(std::is_same_v<typename YVector::value_type, typename YVector::non_const_value_type>,
                "KokkosBlas::swap: YVector must store non const values.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename YVector::memory_space>::accessible,
                "swap: execution_space cannot access data in YVector");

  using XVector_Internal = Kokkos::View<
      typename XVector::non_const_value_type*, typename KokkosKernels::Impl::GetUnifiedLayout<XVector>::array_layout,
      Kokkos::Device<execution_space, typename XVector::memory_space>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
  using YVector_Internal = Kokkos::View<
      typename YVector::non_const_value_type*, typename KokkosKernels::Impl::GetUnifiedLayout<YVector>::array_layout,
      Kokkos::Device<execution_space, typename YVector::memory_space>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  XVector_Internal X(x);
  YVector_Internal Y(y);

  // Runtime check of the length of X and Y
  if (static_cast<int64_t>(X.extent(0)) != static_cast<int64_t>(Y.extent(0))) {
    throw std::runtime_error("X and Y must have equal extents!");
  }

  Kokkos::Profiling::pushRegion("KokkosBlas::swap");
  // If X.extent(0) == 0, do nothing
  if (X.extent(0) != 0) {
    Impl::Swap<execution_space, XVector_Internal, YVector_Internal>::swap(space, X, Y);
  }
  Kokkos::Profiling::popRegion();
}

/// \brief Swaps the entries of vectors x and y.
///
/// \tparam XVector Type of the first vector x; a rank 1 Kokkos::View.
/// \tparam YVector Type of the first vector y; a rank 1 Kokkos::View.
///
/// \param x [in/out] rank 1 View.
/// \param y [in/out] rank 1 View.
///
/// This function is non-blocking unless the underlying TPL requested
/// at compile time is itself blocking. Note that the kernel will be
/// executed on the default stream of the execution_space associted with x.
template <class XVector, class YVector>
void swap(const XVector& x, const YVector& y) {
  const typename XVector::execution_space space = typename XVector::execution_space();
  swap(space, x, y);
}

}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_SWAP_HPP_
