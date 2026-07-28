// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_ROT_HPP_
#define KOKKOSBATCHED_ROT_HPP_

#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Serial Batched Rot:
/// Applies a plane rotation to vectors x and y:
///   x(i) := c*x(i) + s*y(i)
///   y(i) := c*y(i) - s*x(i)          (Conj = false, {s,d,cs,zd}rot)
///   y(i) := c*y(i) - conj(s)*x(i)    (Conj = true, {c,z}rot)
///
/// \tparam Conj: Type indicating whether s is used directly (false)
/// or its conjugate is used (true) in the update of y
///
/// \tparam XViewType: Input/output type for the vector x, needs to be a 1D view
/// \tparam YViewType: Input/output type for the vector y, needs to be a 1D view
/// \tparam CType: Input type for the cosine c (typically real)
/// \tparam SType: Input type for the sine s (real or complex)
///
/// \param[in,out] x: x is a length n vector, a rank 1 view
/// \param[in,out] y: y is a length n vector, a rank 1 view
/// \param[in] c: cosine of the rotation (real scalar)
/// \param[in] s: sine of the rotation (real or complex scalar)
///
/// No nested parallel_for is used inside of the function.
///
template <bool Conj = false>
struct SerialRot {
  template <typename XViewType, typename YViewType, typename CType, typename SType>
  KOKKOS_INLINE_FUNCTION static int invoke(const XViewType &x, const YViewType &y, const CType c, const SType s);
};

/// \brief Team Batched Rot:
/// Applies a plane rotation to vectors x and y:
///   x(i) := c*x(i) + s*y(i)
///   y(i) := c*y(i) - s*x(i)          (Conj = false, {s,d,cs,zd}rot)
///   y(i) := c*y(i) - conj(s)*x(i)    (Conj = true, {c,z}rot)
///
/// \tparam MemberType: TeamPolicy member type
/// \tparam Conj: Type indicating whether s is used directly (false)
/// or its conjugate is used (true) in the update of y
///
/// \tparam XViewType: Input/output type for the vector x, needs to be a 1D view
/// \tparam YViewType: Input/output type for the vector y, needs to be a 1D view
/// \tparam CType: Input type for the cosine c (typically real)
/// \tparam SType: Input type for the sine s (real or complex)
///
/// \param[in] member: TeamPolicy member
/// \param[in,out] x: x is a length n vector, a rank 1 view
/// \param[in,out] y: y is a length n vector, a rank 1 view
/// \param[in] c: cosine of the rotation (real scalar)
/// \param[in] s: sine of the rotation (real or complex scalar)
///
/// A nested parallel_for with TeamThreadRange is used.
///
template <typename MemberType, bool Conj = false>
struct TeamRot {
  template <typename XViewType, typename YViewType, typename CType, typename SType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &x, const YViewType &y,
                                           const CType c, const SType s);
};

/// \brief TeamVector Batched Rot:
/// Applies a plane rotation to vectors x and y:
///   x(i) := c*x(i) + s*y(i)
///   y(i) := c*y(i) - s*x(i)          (Conj = false, {s,d,cs,zd}rot)
///   y(i) := c*y(i) - conj(s)*x(i)    (Conj = true, {c,z}rot)
///
/// \tparam MemberType: TeamPolicy member type
/// \tparam Conj: Type indicating whether s is used directly (false)
/// or its conjugate is used (true) in the update of y
///
/// \tparam XViewType: Input/output type for the vector x, needs to be a 1D view
/// \tparam YViewType: Input/output type for the vector y, needs to be a 1D view
/// \tparam CType: Input type for the cosine c (typically real)
/// \tparam SType: Input type for the sine s (real or complex)
///
/// \param[in] member: TeamPolicy member
/// \param[in,out] x: x is a length n vector, a rank 1 view
/// \param[in,out] y: y is a length n vector, a rank 1 view
/// \param[in] c: cosine of the rotation (real scalar)
/// \param[in] s: sine of the rotation (real or complex scalar)
///
/// A nested parallel_for with TeamVectorRange is used.
///
template <typename MemberType, bool Conj = false>
struct TeamVectorRot {
  template <typename XViewType, typename YViewType, typename CType, typename SType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &x, const YViewType &y,
                                           const CType c, const SType s);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Rot_Impl.hpp"

#endif  // KOKKOSBATCHED_ROT_HPP_
