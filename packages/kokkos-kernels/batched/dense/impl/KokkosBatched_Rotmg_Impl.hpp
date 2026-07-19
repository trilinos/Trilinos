// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_ROTMG_IMPL_HPP_
#define KOKKOSBATCHED_ROTMG_IMPL_HPP_

#include <KokkosBlas_util.hpp>
#include <KokkosBatched_Util.hpp>
#include "KokkosBlas1_rotmg_impl.hpp"

namespace KokkosBatched {
/// \brief invoke Rotmg
/// constructs a modified Givens transformation matrix H which zeros the second component of the 2-vector
/// [sqrt(d1)*x1, sqrt(d2)*y1]**T, where d1, d2, x1, y1 are scalars. The first component of the rotated vector is stored
/// in d1. The elements of H are stored in the array param.
///
/// The matrix H is given by
/// flag == -1.0  flag ==  0.0  flag ==  1.0  flag == -2.0
/// [[h11, h12],  [[1, h12],    [[h11, 1],    [[1, 0],
///  [h21, h22]]   [h21, 1]]     [[-1, h12]]   [[0, 1]]
/// \tparam DXViewType 0-D View type containing a nonconst real scalar
/// \tparam YViewType 0-D View type containing a real scalar
/// \tparam PViewType 1-D View type containing a nonconst real scalar
///
/// \param[in,out] d1: On entry, the scalar d1. On exit, square of x scaling factor to be applied after rotm.
/// \param[in,out] d2: On entry, the scalar d2. On exit, square of y scaling factor to be applied after rotm.
/// \param[in,out] x1: On entry, element from first vector to rotate. On exit, the rotated element before scaling.
/// \param[in] y1: On entry, element from second vector to rotate.
/// \param[out] param: The array containing the elements of the modified Givens transformation matrix H.
///
template <class DXViewType, class YViewType, class PViewType>
KOKKOS_INLINE_FUNCTION int Rotmg::invoke(const DXViewType &d1, const DXViewType &d2, const DXViewType &x1,
                                         const YViewType &y1, const PViewType &param) {
  static_assert(Kokkos::is_view_v<DXViewType>, "KokkosBatched::rotmg: DXViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<YViewType>, "KokkosBatched::rotmg: YViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<PViewType>, "KokkosBatched::rotmg: PViewType is not a Kokkos::View.");
  static_assert(DXViewType::rank() == 0, "KokkosBatched::rotmg: DXViewType must have rank 0.");
  static_assert(YViewType::rank() == 0, "KokkosBatched::rotmg: YViewType must have rank 0.");
  static_assert(PViewType::rank() == 1, "KokkosBatched::rotmg: PViewType must have rank 1.");
  static_assert(std::is_same_v<typename DXViewType::value_type, typename DXViewType::non_const_value_type>,
                "KokkosBatched::rotmg: DXViewType must have non-const value type.");
  static_assert(std::is_same_v<typename YViewType::value_type, typename YViewType::non_const_value_type>,
                "KokkosBatched::rotmg: YViewType must have non-const value type.");
  static_assert(std::is_same_v<typename PViewType::value_type, typename PViewType::non_const_value_type>,
                "KokkosBatched::rotmg: PViewType must have non-const value type.");
  using x_value_type = typename DXViewType::non_const_value_type;
  using y_value_type = typename YViewType::non_const_value_type;
  using p_value_type = typename PViewType::non_const_value_type;
  static_assert(!KokkosKernels::ArithTraits<x_value_type>::is_complex &&
                    !KokkosKernels::ArithTraits<y_value_type>::is_complex &&
                    !KokkosKernels::ArithTraits<p_value_type>::is_complex,
                "KokkosBatched::rotmg: Complex types are not supported for Rotmg.");
#ifndef NDEBUG
  // param should include flag, h11, h21, h12, h22
  if (param.extent_int(0) != 5) {
    Kokkos::printf("KokkosBatched::rotmg: param must have length 5: param length = %d\n", param.extent_int(0));
    return 1;
  }
#endif
  KokkosBlas::Impl::rotmg_impl(d1, d2, x1, y1, param);
  return 0;
}

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_ROTMG_IMPL_HPP_
