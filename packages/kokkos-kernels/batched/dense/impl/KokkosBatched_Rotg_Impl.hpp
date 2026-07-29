// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_ROTG_IMPL_HPP_
#define KOKKOSBATCHED_ROTG_IMPL_HPP_

#include <KokkosBlas_util.hpp>
#include <KokkosBatched_Util.hpp>
#include "KokkosBlas1_rotg_impl.hpp"

namespace KokkosBatched {
/// \brief invoke Rotg for scalars
/// constructs a plane rotation such that
/// [[c,         s],  * [[a],  = [[r],
///  [-conjg(s), c]]     [b]]     [0]]
///
/// \tparam SViewType 0-D View type containing a nonconst real or complex scalar
/// \tparam MViewType 0-D View type containing a nonconst real/magnitude-typed scalar
///
/// \param[in,out] a: On entry, the scalar a. On exit, it is overwritten by r, the nonzero element of the rotated
/// vector.
/// \param[in,out] b: The scalar b (real or complex).
/// \param[out] c: cosine of the rotation (real scalar)
/// \param[out] s: sine of the rotation (real or complex scalar)
///
template <class SViewType, class MViewType>
KOKKOS_INLINE_FUNCTION int Rotg::invoke(const SViewType &a, const SViewType &b, const MViewType &c,
                                        const SViewType &s) {
  static_assert(Kokkos::is_view_v<SViewType>, "KokkosBatched::rotg: SViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<MViewType>, "KokkosBatched::rotg: MViewType is not a Kokkos::View.");
  static_assert(SViewType::rank() == 0, "KokkosBatched::rotg: SViewType must have rank 0.");
  static_assert(MViewType::rank() == 0, "KokkosBatched::rotg: MViewType must have rank 0.");
  static_assert(std::is_same_v<typename SViewType::value_type, typename SViewType::non_const_value_type>,
                "KokkosBatched::rotg: SViewType must have non-const value type.");
  static_assert(std::is_same_v<typename MViewType::value_type, typename MViewType::non_const_value_type>,
                "KokkosBatched::rotg: MViewType must have non-const value type.");
  using MagnitudeType = typename MViewType::non_const_value_type;
  static_assert(!KokkosKernels::ArithTraits<MagnitudeType>::is_complex, "KokkosBatched::rotg: MViewType must be real.");
  KokkosBlas::Impl::rotg_impl(a.data(), b.data(), c.data(), s.data());
  return 0;
}

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_ROTG_IMPL_HPP_
