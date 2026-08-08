// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_ROTG_HPP_
#define KOKKOSBATCHED_ROTG_HPP_

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Device callable Rotg:
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
struct Rotg {
  template <class SViewType, class MViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const SViewType &a, const SViewType &b, const MViewType &c,
                                           const SViewType &s);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Rotg_Impl.hpp"

#endif  // KOKKOSBATCHED_ROTG_HPP_
