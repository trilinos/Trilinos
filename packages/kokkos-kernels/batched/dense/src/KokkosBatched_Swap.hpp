// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_SWAP_HPP_
#define KOKKOSBATCHED_SWAP_HPP_

#include "KokkosBatched_Util.hpp"

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {
/// \brief Serial Batched swap: swap two vectors x and y of length n.
struct SerialSwap {
  /// \tparam XViewType: Type for input X, needs to be a 1D or 2D view
  /// \tparam YViewType: Type for input Y, needs to be a 1D or 2D view
  ///
  /// \param[in,out] x Input/output view of shape (N)
  /// \param[in,out] y Input/output view of shape (N)
  template <typename XViewType, typename YViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const XViewType &x, const YViewType &y);
};

/// \brief Team Batched swap: swap two vectors x and y of length n.
/// A nested parallel_for with TeamThreadRange is used.
/// \tparam MemberType: Kokkos TeamPolicy member type
template <typename MemberType>
struct TeamSwap {
  /// \tparam XViewType: Type for input X, needs to be a 1D or 2D view
  /// \tparam YViewType: Type for input Y, needs to be a 1D or 2D view
  ///
  /// \param[in] member : Kokkos TeamPolicy member type
  /// \param[in,out] x Input/output view of shape (N)
  /// \param[in,out] y Input/output view of shape (N)
  template <typename XViewType, typename YViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &x, const YViewType &y);
};

/// \brief TeamVector Batched swap: swap two vectors x and y of length n.
/// A nested parallel_for with TeamVectorRange is used.
/// \tparam MemberType: Kokkos TeamPolicy member type
template <typename MemberType>
struct TeamVectorSwap {
  /// \tparam XViewType: Type for input X, needs to be a 1D or 2D view
  /// \tparam YViewType: Type for input Y, needs to be a 1D or 2D view
  ///
  /// \param[in] member : Kokkos TeamPolicy member type
  /// \param[in] x Input view of shape (N)
  /// \param[in] y Input view of shape (N)
  template <typename XViewType, typename YViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &x, const YViewType &y);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Swap_Impl.hpp"

#endif  // KOKKOSBATCHED_SWAP_HPP_
