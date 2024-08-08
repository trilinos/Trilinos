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
#ifndef __KOKKOSBATCHED_SVD_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_SVD_SERIAL_IMPL_HPP__

/// \author Brian Kelley (bmkelle@sandia.gov)

#include "KokkosBatched_SVD_Serial_Internal.hpp"

namespace KokkosBatched {
// Version which computes the full factorization
template <typename AViewType, typename UViewType, typename VViewType, typename SViewType, typename WViewType>
KOKKOS_INLINE_FUNCTION int SerialSVD::invoke(SVD_USV_Tag, const AViewType &A, const UViewType &U,
                                             const SViewType &sigma, const VViewType &Vt, const WViewType &work,
                                             typename AViewType::const_value_type tol) {
  static_assert(Kokkos::is_view_v<AViewType> && AViewType::rank == 2, "SVD: A must be a rank-2 view");
  static_assert(Kokkos::is_view_v<UViewType> && UViewType::rank == 2, "SVD: U must be a rank-2 view");
  static_assert(Kokkos::is_view_v<SViewType> && SViewType::rank == 1, "SVD: s must be a rank-1 view");
  static_assert(Kokkos::is_view_v<VViewType> && VViewType::rank == 2, "SVD: V must be a rank-2 view");
  static_assert(Kokkos::is_view_v<WViewType> && WViewType::rank == 1, "SVD: W must be a rank-1 view");
  static_assert(!std::is_same_v<typename WViewType::array_layout, Kokkos::LayoutStride>,
                "SVD: W must be contiguous (not LayoutStride)");
  using value_type = typename AViewType::non_const_value_type;
  return KokkosBatched::SerialSVDInternal::invoke<value_type>(
      A.extent(0), A.extent(1), A.data(), A.stride(0), A.stride(1), U.data(), U.stride(0), U.stride(1), Vt.data(),
      Vt.stride(0), Vt.stride(1), sigma.data(), sigma.stride(0), work.data(), tol);
}

// Version which computes only singular values
template <typename AViewType, typename SViewType, typename WViewType>
KOKKOS_INLINE_FUNCTION int SerialSVD::invoke(SVD_S_Tag, const AViewType &A, const SViewType &sigma,
                                             const WViewType &work, typename AViewType::const_value_type tol) {
  static_assert(Kokkos::is_view_v<AViewType> && AViewType::rank == 2, "SVD: A must be a rank-2 view");
  static_assert(Kokkos::is_view_v<SViewType> && SViewType::rank == 1, "SVD: s must be a rank-1 view");
  static_assert(Kokkos::is_view_v<WViewType> && WViewType::rank == 1, "SVD: W must be a rank-1 view");
  static_assert(!std::is_same_v<typename WViewType::array_layout, Kokkos::LayoutStride>,
                "SVD: W must be contiguous (not LayoutStride)");
  using value_type = typename AViewType::non_const_value_type;
  return KokkosBatched::SerialSVDInternal::invoke<value_type>(A.extent(0), A.extent(1), A.data(), A.stride(0),
                                                              A.stride(1), nullptr, 0, 0, nullptr, 0, 0, sigma.data(),
                                                              sigma.stride(0), work.data(), tol);
}

}  // namespace KokkosBatched

#endif
