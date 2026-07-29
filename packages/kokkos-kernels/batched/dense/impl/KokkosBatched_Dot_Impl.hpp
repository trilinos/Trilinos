// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_DOT_IMPL_HPP
#define KOKKOSBATCHED_DOT_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Dot_Internal.hpp"

namespace KokkosBatched {
namespace Impl {
template <int First, int... Rest>
struct ExtractAxis {
  static constexpr int value = First;
};

template <int Axis, typename XViewType, typename YViewType, typename NormViewType>
KOKKOS_INLINE_FUNCTION static int checkDotInput([[maybe_unused]] const XViewType &x,
                                                [[maybe_unused]] const YViewType &y,
                                                [[maybe_unused]] const NormViewType &dot) {
  static_assert(Kokkos::is_view_v<XViewType>, "KokkosBatched::dot: XViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<YViewType>, "KokkosBatched::dot: YViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<NormViewType>, "KokkosBatched::dot: NormViewType is not a Kokkos::View.");
  static_assert(XViewType::rank() == 1 || XViewType::rank() == 2,
                "KokkosBatched::dot: XViewType must have rank 1 or 2.");
  static_assert(YViewType::rank() == 1 || YViewType::rank() == 2,
                "KokkosBatched::dot: YViewType must have rank 1 or 2.");
  static_assert(NormViewType::rank() == 0 || NormViewType::rank() == 1,
                "KokkosBatched::dot: NormViewType must have rank 0 or 1.");
  static_assert(XViewType::rank() == YViewType::rank() && XViewType::rank() == NormViewType::rank() + 1,
                "KokkosBatched::dot: XViewType and YViewType must have the same rank and must be one rank higher than "
                "NormViewType.");
  static_assert(Axis < XViewType::rank(), "KokkosBatched::dot: Axis must be less than the rank of XViewType.");
#ifndef NDEBUG
  if constexpr (XViewType::rank() == 2) {
    const int x0 = x.extent_int(0), x1 = x.extent_int(1);
    const int y0 = y.extent_int(0), y1 = y.extent_int(1);
    const int dot0 = dot.extent_int(0);

    if (x0 != y0 || x1 != y1) {
      Kokkos::printf(
          "KokkosBatched::dot: Dimensions of X and Y do not match: X: %d x %d, "
          "Y: %d x %d\n",
          x0, x1, y0, y1);
      return 1;
    }

    if (dot0 != (Axis == 0 ? x1 : x0)) {
      Kokkos::printf(
          "KokkosBatched::dot: Dimension of dot does not match the dimension of X and Y: "
          "dot: %d, expected: %d\n",
          dot0, Axis == 0 ? x1 : x0);
      return 1;
    }
  } else {
    if (y.extent_int(0) != x.extent_int(0)) {
      Kokkos::printf(
          "KokkosBatched::dot: x and y must have the same length: x length "
          "= "
          "%d, y length = %d\n",
          x.extent_int(0), y.extent_int(0));
      return 1;
    }
  }
#endif
  return 0;
}
}  // namespace Impl

///
/// Serial Internal Impl
/// ====================
template <typename ArgTrans, int... Args>
template <typename XViewType, typename YViewType, typename NormViewType>
KOKKOS_INLINE_FUNCTION int SerialDot<ArgTrans, Args...>::invoke(const XViewType &X, const YViewType &Y,
                                                                const NormViewType &dot) {
  // Quick return if possible
  const int n = X.size();
  if (n == 0) return 0;

  if constexpr (sizeof...(Args) == 1) {
    constexpr int Axis = Impl::ExtractAxis<Args...>::value;

    auto info = Impl::checkDotInput<Axis>(X, Y, dot);
    if (info) return info;

    using op = std::conditional_t<std::is_same_v<ArgTrans, Trans::ConjTranspose>, KokkosBlas::Impl::OpConj,
                                  KokkosBlas::Impl::OpID>;

    if constexpr (XViewType::rank() == 1) {
      return Impl::SerialDotInternal::invoke(op(), X.extent(0), X.data(), X.stride(0), Y.data(), Y.stride(0),
                                             dot.data());
    } else {
      if constexpr (Axis == 0) {
        // Axis 0 (row-based) : C(j) = op(A(:,j))*B(:,j)
        return Impl::SerialDotInternal::invoke(op(), X.extent(0), X.extent(1), X.data(), X.stride(0), X.stride(1),
                                               Y.data(), Y.stride(0), Y.stride(1), dot.data(), dot.stride(0));

      } else {
        // Axis 1 (column-based) : C(i) = op(A(i,:))*B(i,:)
        return Impl::SerialDotInternal::invoke(op(), X.extent(1), X.extent(0), X.data(), X.stride(1), X.stride(0),
                                               Y.data(), Y.stride(1), Y.stride(0), dot.data(), dot.stride(0));
      }
    }
  } else {
    // Old API
    if constexpr (std::is_same_v<ArgTrans, Trans::Transpose>) {
      // Axis 0 (row-based) : C(j) = conj(A(:,j))*B(:,j)
      return Impl::SerialDotInternal::invoke(KokkosBlas::Impl::OpConj(), X.extent(0), X.extent(1), X.data(),
                                             X.stride(0), X.stride(1), Y.data(), Y.stride(0), Y.stride(1), dot.data(),
                                             dot.stride(0));

    } else {
      // Axis 1 (column-based) : C(i) = conj(A(i,:))*B(i,:)
      return Impl::SerialDotInternal::invoke(KokkosBlas::Impl::OpConj(), X.extent(1), X.extent(0), X.data(),
                                             X.stride(1), X.stride(0), Y.data(), Y.stride(1), Y.stride(0), dot.data(),
                                             dot.stride(0));
    }
  }
}

///
/// Team Impl
/// ===============

template <typename MemberType, typename ArgTrans, int... Args>
template <typename XViewType, typename YViewType, typename NormViewType>
KOKKOS_INLINE_FUNCTION int TeamDot<MemberType, ArgTrans, Args...>::invoke(const MemberType &member, const XViewType &X,
                                                                          const YViewType &Y, const NormViewType &dot) {
  // Quick return if possible
  const int n = X.size();
  if (n == 0) return 0;

  if constexpr (sizeof...(Args) == 1) {
    constexpr int Axis = Impl::ExtractAxis<Args...>::value;
    auto info          = Impl::checkDotInput<Axis>(X, Y, dot);
    if (info) return info;

    using op = std::conditional_t<std::is_same_v<ArgTrans, Trans::ConjTranspose>, KokkosBlas::Impl::OpConj,
                                  KokkosBlas::Impl::OpID>;

    if constexpr (XViewType::rank() == 1) {
      return Impl::TeamDotInternal::invoke(member, op(), X.extent(0), X.data(), X.stride(0), Y.data(), Y.stride(0),
                                           dot.data());
    } else {
      if constexpr (Axis == 0) {
        // Axis 0 (row-based) : C(j) = op(A(:,j))*B(:,j)
        if (X.extent(1) == 1) {
          return Impl::TeamDotInternal::invoke(member, op(), X.extent(0), X.data(), X.stride(0), Y.data(), Y.stride(0),
                                               dot.data());
        }
        return Impl::TeamDotInternal::invoke(member, op(), X.extent(0), X.extent(1), X.data(), X.stride(0), X.stride(1),
                                             Y.data(), Y.stride(0), Y.stride(1), dot.data(), dot.stride(0));
      } else {
        // Axis 1 (column-based) : C(i) = op(A(i,:))*B(i,:)
        if (X.extent(0) == 1) {
          return Impl::TeamDotInternal::invoke(member, op(), X.extent(1), X.data(), X.stride(1), Y.data(), Y.stride(1),
                                               dot.data());
        }
        return Impl::TeamDotInternal::invoke(member, op(), X.extent(1), X.extent(0), X.data(), X.stride(1), X.stride(0),
                                             Y.data(), Y.stride(1), Y.stride(0), dot.data(), dot.stride(0));
      }
    }
  } else {
    // Old API
    if constexpr (std::is_same_v<ArgTrans, Trans::Transpose>) {
      // Axis 0 (row-based) : C(j) = conj(A(:,j))*B(:,j)
      return Impl::TeamDotInternal::invoke(member, KokkosBlas::Impl::OpConj(), X.extent(0), X.extent(1), X.data(),
                                           X.stride(0), X.stride(1), Y.data(), Y.stride(0), Y.stride(1), dot.data(),
                                           dot.stride(0));

    } else {
      // Axis 1 (column-based) : C(i) = conj(A(i,:))*B(i,:)
      return Impl::TeamDotInternal::invoke(member, KokkosBlas::Impl::OpConj(), X.extent(1), X.extent(0), X.data(),
                                           X.stride(1), X.stride(0), Y.data(), Y.stride(1), Y.stride(0), dot.data(),
                                           dot.stride(0));
    }
  }
}
///
/// TeamVector Impl
/// ===============

template <typename MemberType, typename ArgTrans, int... Args>
template <typename XViewType, typename YViewType, typename NormViewType>
KOKKOS_INLINE_FUNCTION int TeamVectorDot<MemberType, ArgTrans, Args...>::invoke(const MemberType &member,
                                                                                const XViewType &X, const YViewType &Y,
                                                                                const NormViewType &dot) {
  // Quick return if possible
  const int n = X.size();
  if (n == 0) return 0;

  if constexpr (sizeof...(Args) == 1) {
    constexpr int Axis = Impl::ExtractAxis<Args...>::value;
    auto info          = Impl::checkDotInput<Axis>(X, Y, dot);
    if (info) return info;

    using op = std::conditional_t<std::is_same_v<ArgTrans, Trans::ConjTranspose>, KokkosBlas::Impl::OpConj,
                                  KokkosBlas::Impl::OpID>;

    if constexpr (XViewType::rank() == 1) {
      return Impl::TeamVectorDotInternal::invoke(member, op(), X.extent(0), X.data(), X.stride(0), Y.data(),
                                                 Y.stride(0), dot.data());
    } else {
      if constexpr (Axis == 0) {
        // Axis 0 (row-based) : C(j) = op(A(:,j))*B(:,j)
        if (X.extent(1) == 1) {
          return Impl::TeamVectorDotInternal::invoke(member, op(), X.extent(0), X.data(), X.stride(0), Y.data(),
                                                     Y.stride(0), dot.data());
        }
        return Impl::TeamVectorDotInternal::invoke(member, op(), X.extent(0), X.extent(1), X.data(), X.stride(0),
                                                   X.stride(1), Y.data(), Y.stride(0), Y.stride(1), dot.data(),
                                                   dot.stride(0));
      } else {
        // Axis 1 (column-based) : C(i) = op(A(i,:))*B(i,:)
        if (X.extent(0) == 1) {
          return Impl::TeamVectorDotInternal::invoke(member, op(), X.extent(1), X.data(), X.stride(1), Y.data(),
                                                     Y.stride(1), dot.data());
        }
        return Impl::TeamVectorDotInternal::invoke(member, op(), X.extent(1), X.extent(0), X.data(), X.stride(1),
                                                   X.stride(0), Y.data(), Y.stride(1), Y.stride(0), dot.data(),
                                                   dot.stride(0));
      }
    }
  } else {
    // Old API
    if constexpr (std::is_same_v<ArgTrans, Trans::Transpose>) {
      // Axis 0 (row-based) : C(j) = conj(A(:,j))*B(:,j)
      return Impl::TeamVectorDotInternal::invoke(member, KokkosBlas::Impl::OpConj(), X.extent(0), X.extent(1), X.data(),
                                                 X.stride(0), X.stride(1), Y.data(), Y.stride(0), Y.stride(1),
                                                 dot.data(), dot.stride(0));

    } else {
      // Axis 1 (column-based) : C(i) = conj(A(i,:))*B(i,:)
      return Impl::TeamVectorDotInternal::invoke(member, KokkosBlas::Impl::OpConj(), X.extent(1), X.extent(0), X.data(),
                                                 X.stride(1), X.stride(0), Y.data(), Y.stride(1), Y.stride(0),
                                                 dot.data(), dot.stride(0));
    }
  }
}

}  // end namespace KokkosBatched

#endif
