// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_DOT_HPP
#define KOKKOSBATCHED_DOT_HPP

/// \author Kim Liegeois (knliege@sandia.gov)
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

template <bool>
struct SerialDot_Deprecated_Warning {};

template <>
struct [[deprecated(
    "KokkosBatched::SerialDot<ArgTrans> is deprecated. Please use "
    "KokkosBatched::SerialDot<ArgTrans, Axis> instead.")]] SerialDot_Deprecated_Warning<true> {};

/// \brief Serial Batched DOT:
///
/// Depending on the Axis template, the dot product is
/// row-based (Axis == 0):
/// C(j) = op(A(:,j))*B(:,j)
///
/// Or column-based (Axis == 1):
/// C(i) = op(A(i,:))*B(i,:)
///
/// For 1D case, C = op(A(:))*B(:)
///
/// No nested parallel_for is used inside of the function.
/// \tparam ArgTrans: whether to apply transpose or conjugate transpose to the input views
/// \tparam Axis: type of dot product (0 for row-based, 1 for column-based)
template <typename ArgTrans = Trans::Transpose, int... Args>
struct SerialDot : SerialDot_Deprecated_Warning<sizeof...(Args) == 0> {
  static_assert(KokkosBlas::is_trans_v<ArgTrans>, "KokkosBatched::SerialDot: ArgTrans must be a KokkosBlas::Trans.");
  /// \tparam XViewType: Input type for X, needs to be a 1D or 2D view
  /// \tparam YViewType: Input type for Y, needs to be a 1D or 2D view
  /// \tparam NormViewType: Input type for alpha, needs to be a 0D or 1D view
  ///
  /// \param[in] X : Input vector X, a rank 1 or 2 view
  /// \param[in] Y : Input vector Y, a rank 1 or 2 view
  /// \param[out] dot : Computed dot product, a rank 0 or 1 view
  template <typename XViewType, typename YViewType, typename NormViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const XViewType &X, const YViewType &Y, const NormViewType &dot);
};

template <bool>
struct TeamDot_Deprecated_Warning {};

template <>
struct [[deprecated(
    "KokkosBatched::TeamDot<ArgTrans> is deprecated. Please use "
    "KokkosBatched::TeamDot<ArgTrans, Axis> instead.")]] TeamDot_Deprecated_Warning<true> {};

/// \brief Team Batched DOT:
///
/// Depending on the Axis template, the dot product is
/// row-based (Axis == 0):
/// C(j) = op(A(:,j))*B(:,j)
///
/// Or column-based (Axis == 1):
/// C(i) = op(A(i,:))*B(i,:)
///
/// For 1D case, C = op(A(:))*B(:)
///
/// \tparam MemberType: Kokkos TeamPolicy member type
/// \tparam ArgTrans: whether to apply transpose or conjugate transpose to the input views
/// \tparam Axis: type of dot product (0 for row-based, 1 for column-based)
template <typename MemberType, typename ArgTrans = Trans::Transpose, int... Args>
struct TeamDot : TeamDot_Deprecated_Warning<sizeof...(Args) == 0> {
  static_assert(KokkosBlas::is_trans_v<ArgTrans>, "KokkosBatched::TeamDot: ArgTrans must be a KokkosBlas::Trans.");
  /// \tparam XViewType: Input type for X, needs to be a 1D or 2D view
  /// \tparam YViewType: Input type for Y, needs to be a 1D or 2D view
  /// \tparam NormViewType: Input type for alpha, needs to be a 0D or 1D view
  ///
  /// \param[in] member : Kokkos TeamPolicy member type
  /// \param[in] X : Input vector X, a rank 1 or 2 view
  /// \param[in] Y : Input vector Y, a rank 1 or 2 view
  /// \param[out] dot : Computed dot product, a rank 0 or 1 view
  ///
  /// A nested parallel_for with TeamThreadRange is used.
  template <typename XViewType, typename YViewType, typename NormViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &X, const YViewType &Y,
                                           const NormViewType &dot);
};

template <bool>
struct TeamVectorDot_Deprecated_Warning {};

template <>
struct [[deprecated(
    "KokkosBatched::TeamVectorDot<ArgTrans> is deprecated. Please use "
    "KokkosBatched::TeamVectorDot<ArgTrans, Axis> instead.")]] TeamVectorDot_Deprecated_Warning<true> {};

/// \brief TeamVector Batched DOT:
///
/// Depending on the Axis template, the dot product is
/// row-based (Axis == 0):
/// C(j) = op(A(:,j))*B(:,j)
///
/// Or column-based (Axis == 1):
/// C(i) = op(A(i,:))*B(i,:)
///
/// For 1D case, C = op(A(:))*B(:)
///
/// Two nested parallel_for with both TeamThreadRange and ThreadVectorRange
/// (or one with TeamVectorRange) are used inside.
///
/// \tparam MemberType: Kokkos TeamPolicy member type
/// \tparam ArgTrans: whether to apply transpose or conjugate transpose to the input views
/// \tparam Axis: type of dot product (0 for row-based, 1 for column-based)
template <typename MemberType, typename ArgTrans = Trans::Transpose, int... Args>
struct TeamVectorDot : TeamVectorDot_Deprecated_Warning<sizeof...(Args) == 0> {
  static_assert(KokkosBlas::is_trans_v<ArgTrans>,
                "KokkosBatched::TeamVectorDot: ArgTrans must be a KokkosBlas::Trans.");
  /// \tparam XViewType: Input type for X, needs to be a 1D or 2D view
  /// \tparam YViewType: Input type for Y, needs to be a 1D or 2D view
  /// \tparam NormViewType: Input type for alpha, needs to be a 0D or 1D view
  ///
  /// \param[in] X : Input vector X, a rank 1 or 2 view
  /// \param[in] Y : Input vector Y, a rank 1 or 2 view
  /// \param[out] dot : Computed dot product, a rank 0 or 1 view
  template <typename XViewType, typename YViewType, typename NormViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &X, const YViewType &Y,
                                           const NormViewType &dot);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Dot_Impl.hpp"

#endif
