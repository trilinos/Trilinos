// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_IAMAX_INTERNAL_HPP_
#define KOKKOSBATCHED_IAMAX_INTERNAL_HPP_

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include <Kokkos_Core.hpp>
#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {
namespace Impl {

///
/// Serial Internal Impl
/// ========================

struct SerialIamaxInternal {
  template <typename IndexType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static IndexType invoke(const IndexType n, const ValueType *KOKKOS_RESTRICT x,
                                                 const IndexType xs0);
};

template <typename IndexType, typename ValueType>
KOKKOS_INLINE_FUNCTION IndexType SerialIamaxInternal::invoke(const IndexType n, const ValueType *KOKKOS_RESTRICT x,
                                                             const IndexType xs0) {
  // Quick return
  if (n <= 1) return 0;

  using ats      = typename KokkosKernels::ArithTraits<ValueType>;
  using RealType = typename ats::mag_type;

  RealType amax  = Kokkos::abs(x[0 * xs0]);
  IndexType imax = 0;

  for (IndexType i = 1; i < n; ++i) {
    const RealType abs_x_i = Kokkos::abs(x[i * xs0]);
    if (abs_x_i > amax) {
      amax = abs_x_i;
      imax = i;
    }
  }

  return imax;
}

///
/// Team Internal Impl
/// ========================
template <typename MemberType>
struct TeamIamaxInternal {
  template <typename IndexType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static IndexType invoke(const MemberType &member, const IndexType n,
                                                 const ValueType *KOKKOS_RESTRICT x, const IndexType xs0);
};

template <typename MemberType>
template <typename IndexType, typename ValueType>
KOKKOS_INLINE_FUNCTION IndexType TeamIamaxInternal<MemberType>::invoke(const MemberType &member, const IndexType n,
                                                                       const ValueType *KOKKOS_RESTRICT x,
                                                                       const IndexType xs0) {
  // Quick return
  if (n <= 1) return 0;

  using mag_type           = typename KokkosKernels::ArithTraits<ValueType>::mag_type;
  using reducer_type       = Kokkos::MaxFirstLoc<mag_type, IndexType, typename MemberType::execution_space>;
  using reducer_value_type = typename reducer_type::value_type;
  reducer_value_type max_first_loc{};
  Kokkos::parallel_reduce(
      Kokkos::TeamThreadRange(member, n),
      [&](const IndexType i, reducer_value_type &update) {
        const auto abs_val = KokkosKernels::ArithTraits<ValueType>::abs(x[i * xs0]);
        if (abs_val > update.val) {
          update.val = abs_val;
          update.loc = i;
        }
      },
      reducer_type(max_first_loc));
  return max_first_loc.loc;
}

///
/// TeamVector Internal Impl
/// ========================
template <typename MemberType>
struct TeamVectorIamaxInternal {
  template <typename IndexType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static IndexType invoke(const MemberType &member, const IndexType n,
                                                 const ValueType *KOKKOS_RESTRICT x, const IndexType xs0);
};

template <typename MemberType>
template <typename IndexType, typename ValueType>
KOKKOS_INLINE_FUNCTION IndexType TeamVectorIamaxInternal<MemberType>::invoke(const MemberType &member,
                                                                             const IndexType n,
                                                                             const ValueType *KOKKOS_RESTRICT x,
                                                                             const IndexType xs0) {
  // Quick return
  if (n <= 1) return 0;

  using mag_type           = typename KokkosKernels::ArithTraits<ValueType>::mag_type;
  using reducer_type       = Kokkos::MaxFirstLoc<mag_type, IndexType, typename MemberType::execution_space>;
  using reducer_value_type = typename reducer_type::value_type;
  reducer_value_type max_first_loc{};
  Kokkos::parallel_reduce(
      Kokkos::TeamVectorRange(member, n),
      [&](const IndexType i, reducer_value_type &update) {
        const auto abs_val = KokkosKernels::ArithTraits<ValueType>::abs(x[i * xs0]);
        if (abs_val > update.val) {
          update.val = abs_val;
          update.loc = i;
        }
      },
      reducer_type(max_first_loc));
  return max_first_loc.loc;
}

}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_IAMAX_INTERNAL_HPP_
