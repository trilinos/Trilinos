// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_COPY_IMPL_HPP
#define KOKKOSBATCHED_COPY_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include "KokkosBlas_util.hpp"
#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Copy_Internal.hpp"

namespace KokkosBatched {
namespace Impl {

template <typename ArgTrans, typename AViewType, typename BViewType>
KOKKOS_INLINE_FUNCTION static int checkCopyInput([[maybe_unused]] const AViewType &A,
                                                 [[maybe_unused]] const BViewType &B) {
  static_assert(Kokkos::is_view_v<AViewType>, "KokkosBatched::copy: AViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<BViewType>, "KokkosBatched::copy: BViewType is not a Kokkos::View.");
  static_assert(AViewType::rank() == BViewType::rank(),
                "KokkosBatched::copy: AViewType and BViewType must have the same rank.");
  static_assert(AViewType::rank() == 1 || AViewType::rank() == 2,
                "KokkosBatched::copy: Only rank 1 and rank 2 views are supported.");
  static_assert(std::is_same_v<typename BViewType::value_type, typename BViewType::non_const_value_type>,
                "KokkosBatched::copy: The value_type of BViewType must have non-const value type.");

#ifndef NDEBUG
  if constexpr (AViewType::rank() == 1) {
    if (A.extent_int(0) != B.extent_int(0)) {
      Kokkos::printf(
          "KokkosBatched::copy: Dimensions of A and B do not match: A: %d, "
          "B: %d\n",
          A.extent_int(0), B.extent_int(0));
      return 1;
    }
  } else {
    const int m_A = std::is_same_v<ArgTrans, Trans::NoTranspose> ? A.extent_int(0) : A.extent_int(1);
    const int n_A = std::is_same_v<ArgTrans, Trans::NoTranspose> ? A.extent_int(1) : A.extent_int(0);
    if (m_A != B.extent_int(0) || n_A != B.extent_int(1)) {
      Kokkos::printf(
          "KokkosBatched::copy: Dimensions of A and B do not match: A: %d x %d, "
          "B: %d x %d\n",
          A.extent_int(0), A.extent_int(1), B.extent_int(0), B.extent_int(1));
      return 1;
    }
  }
#endif
  return 0;
}
}  // namespace Impl

///
/// Serial Impl
/// ===========

template <>
template <typename AViewType, typename BViewType>
KOKKOS_INLINE_FUNCTION int SerialCopy<Trans::NoTranspose>::invoke(const AViewType &A, const BViewType &B) {
  auto info = Impl::checkCopyInput<Trans::NoTranspose>(A, B);
  if (info) return info;

  // Quick return if possible
  if (A.size() == 0 || B.size() == 0) return 0;

  if constexpr (AViewType::rank() == 1) {
    return Impl::SerialCopyInternal::invoke(KokkosBlas::Impl::OpID(), A.extent(0), A.data(), A.stride(0), B.data(),
                                            B.stride(0));
  } else {
    return Impl::SerialCopyInternal::invoke(KokkosBlas::Impl::OpID(), A.extent(0), A.extent(1), A.data(), A.stride(0),
                                            A.stride(1), B.data(), B.stride(0), B.stride(1));
  }
}

template <>
template <typename AViewType, typename BViewType>
KOKKOS_INLINE_FUNCTION int SerialCopy<Trans::Transpose>::invoke(const AViewType &A, const BViewType &B) {
  auto info = Impl::checkCopyInput<Trans::Transpose>(A, B);
  if (info) return info;

  // Quick return if possible
  if (A.size() == 0 || B.size() == 0) return 0;

  if constexpr (AViewType::rank() == 1) {
    return Impl::SerialCopyInternal::invoke(KokkosBlas::Impl::OpID(), A.extent(0), A.data(), A.stride(0), B.data(),
                                            B.stride(0));
  } else {
    return Impl::SerialCopyInternal::invoke(KokkosBlas::Impl::OpID(), A.extent(1), A.extent(0), A.data(), A.stride(1),
                                            A.stride(0), B.data(), B.stride(0), B.stride(1));
  }
}

template <>
template <typename AViewType, typename BViewType>
KOKKOS_INLINE_FUNCTION int SerialCopy<Trans::ConjTranspose>::invoke(const AViewType &A, const BViewType &B) {
  auto info = Impl::checkCopyInput<Trans::ConjTranspose>(A, B);
  if (info) return info;

  // Quick return if possible
  if (A.size() == 0 || B.size() == 0) return 0;

  if constexpr (AViewType::rank() == 1) {
    return Impl::SerialCopyInternal::invoke(KokkosBlas::Impl::OpConj(), A.extent(0), A.data(), A.stride(0), B.data(),
                                            B.stride(0));
  } else {
    return Impl::SerialCopyInternal::invoke(KokkosBlas::Impl::OpConj(), A.extent(1), A.extent(0), A.data(), A.stride(1),
                                            A.stride(0), B.data(), B.stride(0), B.stride(1));
  }
}

///
/// Team Impl
/// =========
template <typename MemberType>
struct TeamCopy<MemberType, Trans::NoTranspose> {
  template <typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B) {
    auto info = Impl::checkCopyInput<Trans::NoTranspose>(A, B);
    if (info) return info;

    // Quick return if possible
    if (A.size() == 0 || B.size() == 0) return 0;

    if constexpr (AViewType::rank() == 1) {
      return Impl::TeamCopyInternal::invoke(member, KokkosBlas::Impl::OpID(), A.extent(0), A.data(), A.stride(0),
                                            B.data(), B.stride(0));
    } else {
      if (A.extent(0) == 1) {
        return Impl::TeamCopyInternal::invoke(member, KokkosBlas::Impl::OpID(), A.extent(1), A.data(), A.stride(1),
                                              B.data(), B.stride(1));
      }
      return Impl::TeamCopyInternal::invoke(member, KokkosBlas::Impl::OpID(), A.extent(0), A.extent(1), A.data(),
                                            A.stride(0), A.stride(1), B.data(), B.stride(0), B.stride(1));
    }
  }
};

template <typename MemberType>
struct TeamCopy<MemberType, Trans::Transpose> {
  template <typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B) {
    auto info = Impl::checkCopyInput<Trans::Transpose>(A, B);
    if (info) return info;

    // Quick return if possible
    if (A.size() == 0 || B.size() == 0) return 0;

    if constexpr (AViewType::rank() == 1) {
      return Impl::TeamCopyInternal::invoke(member, KokkosBlas::Impl::OpID(), A.extent(0), A.data(), A.stride(0),
                                            B.data(), B.stride(0));
    } else {
      if (A.extent(1) == 1) {
        return Impl::TeamCopyInternal::invoke(member, KokkosBlas::Impl::OpID(), A.extent(0), A.data(), A.stride(0),
                                              B.data(), B.stride(1));
      }
      return Impl::TeamCopyInternal::invoke(member, KokkosBlas::Impl::OpID(), A.extent(1), A.extent(0), A.data(),
                                            A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1));
    }
  }
};

template <typename MemberType>
struct TeamCopy<MemberType, Trans::ConjTranspose> {
  template <typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B) {
    auto info = Impl::checkCopyInput<Trans::ConjTranspose>(A, B);
    if (info) return info;

    // Quick return if possible
    if (A.size() == 0 || B.size() == 0) return 0;

    if constexpr (AViewType::rank() == 1) {
      return Impl::TeamCopyInternal::invoke(member, KokkosBlas::Impl::OpConj(), A.extent(0), A.data(), A.stride(0),
                                            B.data(), B.stride(0));
    } else {
      if (A.extent(1) == 1) {
        return Impl::TeamCopyInternal::invoke(member, KokkosBlas::Impl::OpConj(), A.extent(0), A.data(), A.stride(0),
                                              B.data(), B.stride(1));
      }
      return Impl::TeamCopyInternal::invoke(member, KokkosBlas::Impl::OpConj(), A.extent(1), A.extent(0), A.data(),
                                            A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1));
    }
  }
};

///
/// TeamVector Impl
/// =========

template <typename MemberType>
struct TeamVectorCopy<MemberType, Trans::NoTranspose> {
  template <typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B) {
    auto info = Impl::checkCopyInput<Trans::NoTranspose>(A, B);
    if (info) return info;

    // Quick return if possible
    if (A.size() == 0 || B.size() == 0) return 0;

    if constexpr (AViewType::rank() == 1) {
      return Impl::TeamVectorCopyInternal::invoke(member, KokkosBlas::Impl::OpID(), A.extent(0), A.data(), A.stride(0),
                                                  B.data(), B.stride(0));
    } else {
      if (A.extent(0) == 1) {
        return Impl::TeamVectorCopyInternal::invoke(member, KokkosBlas::Impl::OpID(), A.extent(1), A.data(),
                                                    A.stride(1), B.data(), B.stride(1));
      }
      return Impl::TeamVectorCopyInternal::invoke(member, KokkosBlas::Impl::OpID(), A.extent(0), A.extent(1), A.data(),
                                                  A.stride(0), A.stride(1), B.data(), B.stride(0), B.stride(1));
    }
  }
};

template <typename MemberType>
struct TeamVectorCopy<MemberType, Trans::Transpose> {
  template <typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B) {
    auto info = Impl::checkCopyInput<Trans::Transpose>(A, B);
    if (info) return info;

    // Quick return if possible
    if (A.size() == 0 || B.size() == 0) return 0;

    if constexpr (AViewType::rank() == 1) {
      return Impl::TeamVectorCopyInternal::invoke(member, KokkosBlas::Impl::OpID(), A.extent(0), A.data(), A.stride(0),
                                                  B.data(), B.stride(0));
    } else {
      if (A.extent(1) == 1) {
        return Impl::TeamVectorCopyInternal::invoke(member, KokkosBlas::Impl::OpID(), A.extent(0), A.data(),
                                                    A.stride(0), B.data(), B.stride(1));
      }
      return Impl::TeamVectorCopyInternal::invoke(member, KokkosBlas::Impl::OpID(), A.extent(1), A.extent(0), A.data(),
                                                  A.stride(1), A.stride(0), B.data(), B.stride(0), B.stride(1));
    }
  }
};

template <typename MemberType>
struct TeamVectorCopy<MemberType, Trans::ConjTranspose> {
  template <typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const BViewType &B) {
    auto info = Impl::checkCopyInput<Trans::ConjTranspose>(A, B);
    if (info) return info;

    // Quick return if possible
    if (A.size() == 0 || B.size() == 0) return 0;

    if constexpr (AViewType::rank() == 1) {
      return Impl::TeamVectorCopyInternal::invoke(member, KokkosBlas::Impl::OpConj(), A.extent(0), A.data(),
                                                  A.stride(0), B.data(), B.stride(0));
    } else {
      if (A.extent(1) == 1) {
        return Impl::TeamVectorCopyInternal::invoke(member, KokkosBlas::Impl::OpConj(), A.extent(0), A.data(),
                                                    A.stride(0), B.data(), B.stride(1));
      }
      return Impl::TeamVectorCopyInternal::invoke(member, KokkosBlas::Impl::OpConj(), A.extent(1), A.extent(0),
                                                  A.data(), A.stride(1), A.stride(0), B.data(), B.stride(0),
                                                  B.stride(1));
    }
  }
};

}  // end namespace KokkosBatched

#endif
