// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_IAMAX_HPP_
#define KOKKOSBATCHED_IAMAX_HPP_

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Serial Batched Iamax:
/// Iamax finds the index of the first element having maximum absolute value.
///
struct SerialIamax {
  /// \tparam XViewType: Input view type, needs to be a 1D view
  ///
  /// \param[in] X: Input view type
  ///
  /// \return The index of the first element having maximum absolute value
  /// As well as Blas, this returns 0 for an empty vector
  /// No nested parallel_for is used inside of the function.
  ///
  template <typename XViewType>
  KOKKOS_INLINE_FUNCTION static typename XViewType::size_type invoke(const XViewType &x);
};

/// \brief Team Batched Iamax:
/// Iamax finds the index of the first element having maximum absolute value.
///
/// \tparam MemberType: Kokkos TeamPolicy member type
///
template <typename MemberType>
struct TeamIamax {
  /// \tparam XViewType: Input view type, needs to be a 1D view
  /// \param[in] member : Kokkos TeamPolicy member type
  /// \param[in] X: Input view type
  /// \return The index of the first element having maximum absolute value
  /// As well as Blas, this returns 0 for an empty vector
  /// A nested parallel_for with TeamThreadRange is used.
  template <typename XViewType>
  KOKKOS_INLINE_FUNCTION static typename XViewType::size_type invoke(const MemberType &member, const XViewType &x);
};

/// \brief Team Vector Batched Iamax:
/// Iamax finds the index of the first element having maximum absolute value.
///
/// \tparam MemberType: Kokkos TeamPolicy member type
///
template <typename MemberType>
struct TeamVectorIamax {
  /// \tparam XViewType: Input view type, needs to be a 1D view
  /// \param[in] member : Kokkos TeamPolicy member type
  /// \param[in] X: Input view type
  /// \return The index of the first element having maximum absolute value
  /// As well as Blas, this returns 0 for an empty vector
  /// A nested parallel_for with TeamThreadRange is used.
  template <typename XViewType>
  KOKKOS_INLINE_FUNCTION static typename XViewType::size_type invoke(const MemberType &member, const XViewType &x);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Iamax_Impl.hpp"

#endif  // KOKKOSBATCHED_IAMAX_HPP_
