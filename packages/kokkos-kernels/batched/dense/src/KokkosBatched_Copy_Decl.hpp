// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_COPY_DECL_HPP
#define KOKKOSBATCHED_COPY_DECL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

template <bool>
struct SerialCopy_Deprecated_Warning {};

template <>
struct [[deprecated(
    "KokkosBatched::SerialCopy<ArgTrans, rank> is deprecated. Please use "
    "KokkosBatched::SerialCopy<ArgTrans> instead.")]] SerialCopy_Deprecated_Warning<true> {};

/// \brief Serial Batched Copy:
/// Performs B = Op(A)
/// where Op is one of NoTranspose, Transpose, ConjTranspose
/// A and B are 1D or 2D views
///
/// \tparam ArgTrans: one of NoTranspose, Transpose, ConjTranspose
/// \tparam Args: (deprecated) rank information
template <typename ArgTrans = Trans::NoTranspose, int... Args>
struct SerialCopy : SerialCopy_Deprecated_Warning<sizeof...(Args) != 0> {
  static_assert(KokkosBlas::is_trans_v<ArgTrans>, "KokkosBatched::SerialCopy: ArgTrans must be a KokkosBlas::Trans.");

  /// \brief invoke the SerialCopy
  /// \tparam AViewType: Kokkos::View type for A
  /// \tparam BViewType: Kokkos::View type for B
  /// \param[in] A Input view A
  /// \param[out] B Output view B
  template <typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const BViewType &B);
};

template <bool>
struct TeamCopy_Deprecated_Warning {};

template <>
struct [[deprecated(
    "KokkosBatched::TeamCopy<ArgTrans, rank> is deprecated. Please use "
    "KokkosBatched::TeamCopy<ArgTrans> instead.")]] TeamCopy_Deprecated_Warning<true> {};

/// \brief Team Batched Copy:
/// Performs B = Op(A)
/// where Op is one of NoTranspose, Transpose, ConjTranspose
/// A and B are 1D or 2D views
///
/// \tparam MemberType: Kokkos::TeamPolicy member type
/// \tparam ArgTrans: one of NoTranspose, Transpose, ConjTranspose
/// \tparam Args: (deprecated) rank information
template <typename MemberType, typename ArgTrans = Trans::NoTranspose, int... Args>
struct TeamCopy : TeamCopy_Deprecated_Warning<sizeof...(Args) != 0> {
  static_assert(KokkosBlas::is_trans_v<ArgTrans>, "KokkosBatched::TeamCopy: ArgTrans must be a KokkosBlas::Trans.");

  /// \brief invoke the TeamCopy
  /// \tparam AViewType: Kokkos::View type for A
  /// \tparam BViewType: Kokkos::View type for B
  /// \param[in] member Kokkos::TeamPolicy member
  /// \param[in] A Input view A
  /// \param[out] B Output view B
  template <typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B);
};

template <bool>
struct TeamVectorCopy_Deprecated_Warning {};

template <>
struct [[deprecated(
    "KokkosBatched::TeamVectorCopy<ArgTrans, rank> is deprecated. Please use "
    "KokkosBatched::TeamVectorCopy<ArgTrans> instead.")]] TeamVectorCopy_Deprecated_Warning<true> {};

/// \brief TeamVector Batched Copy:
/// Performs B = Op(A)
/// where Op is one of NoTranspose, Transpose, ConjTranspose
/// A and B are 1D or 2D views
///
/// \tparam MemberType: Kokkos::TeamPolicy member type
/// \tparam ArgTrans: one of NoTranspose, Transpose, ConjTranspose
/// \tparam Args: (deprecated) rank information
template <typename MemberType, typename ArgTrans = Trans::NoTranspose, int... Args>
struct TeamVectorCopy : TeamVectorCopy_Deprecated_Warning<sizeof...(Args) != 0> {
  static_assert(KokkosBlas::is_trans_v<ArgTrans>,
                "KokkosBatched::TeamVectorCopy: ArgTrans must be a KokkosBlas::Trans.");

  /// \brief invoke the TeamVectorCopy
  /// \tparam AViewType: Kokkos::View type for A
  /// \tparam BViewType: Kokkos::View type for B
  /// \param[in] member Kokkos::TeamPolicy member
  /// \param[in] A Input view A
  /// \param[out] B Output view B
  template <typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B);
};

template <bool>
struct Copy_Deprecated_Warning {};

template <>
struct [[deprecated(
    "KokkosBatched::Copy<MemberType, ArgTrans, rank> is deprecated. Please use "
    "KokkosBatched::Copy<MemberType, ArgTrans> instead.")]] Copy_Deprecated_Warning<true> {};

/// \brief General Copy:
/// Performs B = Op(A)
/// where Op is one of NoTranspose, Transpose, ConjTranspose
/// A and B are 1D or 2D views
/// Dispatches to SerialCopy, TeamCopy, or TeamVectorCopy based on ArgMode
/// \tparam MemberType: Kokkos::TeamPolicy member type
/// \tparam ArgTrans: one of NoTranspose, Transpose, ConjTranspose
/// \tparam ArgMode: one of Mode::Serial, Mode::Team, Mode::TeamVector
/// \tparam Args: (deprecated) rank information
template <typename MemberType, typename ArgTrans, typename ArgMode, int... Args>
struct Copy : Copy_Deprecated_Warning<sizeof...(Args) != 0> {
  static_assert(KokkosBlas::is_trans_v<ArgTrans>, "KokkosBatched::Copy: ArgTrans must be a KokkosBlas::Trans.");

  /// \brief invoke the Copy
  /// \tparam AViewType: Kokkos::View type for A
  /// \tparam BViewType: Kokkos::View type for B
  /// \param[in] member Kokkos::TeamPolicy member
  /// \param[in] A Input view A
  /// \param[out] B Output view B
  template <typename AViewType, typename BViewType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B) {
    int r_val = 0;
    if constexpr (std::is_same_v<ArgMode, Mode::Serial>) {
      r_val = SerialCopy<ArgTrans>::invoke(A, B);
    } else if constexpr (std::is_same_v<ArgMode, Mode::Team>) {
      r_val = TeamCopy<MemberType, ArgTrans>::invoke(member, A, B);
    } else if constexpr (std::is_same_v<ArgMode, Mode::TeamVector>) {
      r_val = TeamVectorCopy<MemberType, ArgTrans>::invoke(member, A, B);
    }
    return r_val;
  }
};

}  // namespace KokkosBatched

#include "KokkosBatched_Copy_Impl.hpp"

#define KOKKOSBATCHED_SERIAL_COPY_MATRIX_NO_TRANSPOSE_INTERNAL_INVOKE(M, N, A, AS0, AS1, B, BS0, BS1) \
  KokkosBatched::SerialCopyInternal ::invoke(M, N, A, AS0, AS1, B, BS0, BS1)

#define KOKKOSBATCHED_TEAM_COPY_MATRIX_NO_TRANSPOSE_INTERNAL_INVOKE(MEMBER, M, N, A, AS0, AS1, B, BS0, BS1) \
  KokkosBatched::TeamCopyInternal ::invoke(MEMBER, M, N, A, AS0, AS1, B, BS0, BS1)

#define KOKKOSBATCHED_SERIAL_COPY_VECTOR_INTERNAL_INVOKE(M, A, AS, B, BS) \
  KokkosBatched::SerialCopyInternal ::invoke(M, A, AS, B, BS)

#define KOKKOSBATCHED_TEAM_COPY_VECTOR_NO_TRANSPOSE_INTERNAL_INVOKE(MEMBER, M, A, AS, B, BS) \
  KokkosBatched::TeamCopyInternal ::invoke(MEMBER, M, A, AS, B, BS)

#define KOKKOSBATCHED_COPY_VECTOR_NO_TRANSPOSE_INTERNAL_INVOKE(MODETYPE, MEMBER, M, A, AS, B, BS) \
  if constexpr (std::is_same_v<MODETYPE, KokkosBatched::Mode::Serial>) {                          \
    KOKKOSBATCHED_SERIAL_COPY_VECTOR_INTERNAL_INVOKE(M, A, AS, B, BS);                            \
  } else if constexpr (std::is_same_v<MODETYPE, KokkosBatched::Mode::Team>) {                     \
    KOKKOSBATCHED_TEAM_COPY_VECTOR_NO_TRANSPOSE_INTERNAL_INVOKE(MEMBER, M, A, AS, B, BS);         \
  }

#endif
