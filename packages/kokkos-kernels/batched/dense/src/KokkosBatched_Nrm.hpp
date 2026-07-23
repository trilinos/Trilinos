// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_NRM_HPP_
#define KOKKOSBATCHED_NRM_HPP_

#include "KokkosBatched_Util.hpp"

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {
/// \brief Serial Batched nrm:
/// If NrmType == Norm::L1, compute L1 norm of each vector in the batch
/// If NrmType == Norm::L2, compute L2 norm of each vector in the batch
/// If NrmType == Norm::LInf, compute Linf norm of each vector in the batch
/// If NrmType == Norm::ScaledL2, compute ScaledL2 norm of each vector in the batch
/// No nested parallel_for is used inside of the function.
/// \tparam NrmType: one of Norm::L1, Norm::L2, Norm::LInf, Norm::ScaledL2
template <typename NrmType>
struct SerialNrm {
  static_assert(is_norm_v<NrmType>,
                "KokkosBatched::SerialNrm: NrmType must be one of Norm::L1, Norm::L2, Norm::LInf, Norm::ScaledL2");
  /// \tparam XViewType: Type for input X, needs to be a 1D view
  /// \tparam NormViewType: Type for output norm, needs to be a 0D view
  ///
  /// \param[in] x Input view of shape (N)
  /// \param[out] norm Output view of shape () to store the computed norms
  template <typename XViewType, typename NormViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const XViewType &x, const NormViewType &norm);
};

/// \brief Team Batched nrm:
/// If NrmType == Norm::L1, compute L1 norm of each vector in the batch
/// If NrmType == Norm::L2, compute L2 norm of each vector in the batch
/// If NrmType == Norm::LInf, compute Linf norm of each vector in the batch
/// If NrmType == Norm::ScaledL2, compute ScaledL2 norm of each vector in the batch
/// A nested parallel_for with TeamThreadRange is used.
/// \tparam MemberType: Kokkos TeamPolicy member type
/// \tparam NrmType: one of Norm::L1, Norm::L2, Norm::LInf
template <typename MemberType, typename NrmType>
struct TeamNrm {
  static_assert(is_norm_v<NrmType>,
                "KokkosBatched::TeamNrm: NrmType must be one of Norm::L1, Norm::L2, Norm::LInf, Norm::ScaledL2");
  /// \tparam XViewType: Type for input X, needs to be a 1D view
  /// \tparam NormViewType: Type for output norm, needs to be a 0D view
  ///
  /// \param[in] member : Kokkos TeamPolicy member type
  /// \param[in] x Input view of shape (N)
  /// \param[out] norm Output view of shape () to store the computed norms
  template <typename XViewType, typename NormViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &x, const NormViewType &norm);
};

/// \brief Team Batched nrm:
/// If NrmType == Norm::L1, compute L1 norm of each vector in the batch
/// If NrmType == Norm::L2, compute L2 norm of each vector in the batch
/// If NrmType == Norm::LInf, compute Linf norm of each vector in the batch
/// If NrmType == Norm::ScaledL2, compute ScaledL2 norm of each vector in the batch
/// A nested parallel_for with TeamVectorRange is used.
/// \tparam MemberType: Kokkos TeamPolicy member type
/// \tparam NrmType: one of Norm::L1, Norm::L2, Norm::LInf, Norm::ScaledL2
template <typename MemberType, typename NrmType>
struct TeamVectorNrm {
  static_assert(is_norm_v<NrmType>,
                "KokkosBatched::TeamVectorNrm: NrmType must be one of Norm::L1, Norm::L2, Norm::LInf, Norm::ScaledL2");
  /// \tparam XViewType: Type for input X, needs to be a 1D view
  /// \tparam NormViewType: Type for output norm, needs to be a 0D view
  ///
  /// \param[in] member : Kokkos TeamPolicy member type
  /// \param[in] x Input view of shape (N)
  /// \param[out] norm Output view of shape () to store the computed norms
  template <typename XViewType, typename NormViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &x, const NormViewType &norm);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Nrm_Impl.hpp"

#endif  // KOKKOSBATCHED_NRM_HPP_
