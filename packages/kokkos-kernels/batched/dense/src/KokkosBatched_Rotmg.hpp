// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_ROTMG_HPP_
#define KOKKOSBATCHED_ROTMG_HPP_

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Device callable Rotmg:
/// constructs a modified Givens transformation matrix H which zeros the second component of the 2-vector
/// [sqrt(d1)*x1, sqrt(d2)*y1]**T, where d1, d2, x1, y1 are scalars. The first component of the rotated vector is stored
/// in d1. The elements of H are stored in the array param.
///
/// The matrix H is given by
/// flag == -1.0  flag ==  0.0  flag ==  1.0  flag == -2.0
/// [[h11, h12],  [[1, h12],    [[h11, 1],    [[1, 0],
///  [h21, h22]]   [h21, 1]]     [[-1, h12]]   [[0, 1]]
struct Rotmg {
  /// \tparam DXViewType 0-D View type containing a nonconst real scalar
  /// \tparam YViewType 0-D View type containing a real scalar
  /// \tparam PViewType 1-D View type containing a nonconst real scalar
  ///
  /// \param[in,out] d1: On entry, the scalar d1. On exit, square of x scaling factor to be applied after rotm.
  /// \param[in,out] d2: On entry, the scalar d2. On exit, square of y scaling factor to be applied after rotm.
  /// \param[in,out] x1: On entry, element from first vector to rotate. On exit, the rotated element before scaling.
  /// \param[in] y1: On entry, element from second vector to rotate.
  /// \param[out] param: The array containing the elements of the modified Givens transformation matrix H.
  template <class DXViewType, class YViewType, class PViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const DXViewType &d1, const DXViewType &d2, const DXViewType &x1,
                                           const YViewType &y1, const PViewType &param);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Rotmg_Impl.hpp"

#endif  // KOKKOSBATCHED_ROTMG_HPP_
