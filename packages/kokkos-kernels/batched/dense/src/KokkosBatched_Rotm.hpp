// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_ROTM_HPP_
#define KOKKOSBATCHED_ROTM_HPP_

#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Serial Batched Rotm:
/// Applies the modified Givens rotation to vectors x and y:
///  x(i) := h11*x(i) + h12*y(i)
///  y(i) := h21*x(i) + h22*y(i)
///
/// The matrix H is given by
/// flag == -1.0  flag ==  0.0  flag ==  1.0  flag == -2.0
/// [[h11, h12],  [[1, h12],    [[h11, 1],    [[1, 0],
///  [h21, h22]]   [h21, 1]]     [[-1, h12]]   [[0, 1]]
///
/// param is a length 5 vector containing the parameters of the modified Givens rotation:
/// param(0) = flag, param(1) = h11, param(2) = h21, param(3) = h12, param(4) = h22
///
struct SerialRotm {
  /// \brief Invokes the SerialRotm functor. No nested parallel_for is used inside of the function.
  /// \tparam XViewType: Input/output type for the vector x, needs to be a 1D view
  /// \tparam YViewType: Input/output type for the vector y, needs to be a 1D view
  /// \tparam ParamViewType: Input type for the param vector, needs to be a 1D view of length 5
  ///
  /// \param[in,out] x: x is a length n vector, a rank 1 view
  /// \param[in,out] y: y is a length n vector, a rank 1 view
  /// \param[in] param: param is a length 5 vector, a rank 1 view, containing the parameters of the modified Givens
  /// rotation
  template <typename XViewType, typename YViewType, typename ParamViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const XViewType &x, const YViewType &y, const ParamViewType &param);
};

/// \brief Team Batched Rotm:
/// Applies the modified Givens rotation to vectors x and y:
///  x(i) := h11*x(i) + h12*y(i)
///  y(i) := h21*x(i) + h22*y(i)
///
/// The matrix H is given by
/// flag == -1.0  flag ==  0.0  flag ==  1.0  flag == -2.0
/// [[h11, h12],  [[1, h12],    [[h11, 1],    [[1, 0],
///  [h21, h22]]   [h21, 1]]     [[-1, h12]]   [[0, 1]]
///
/// param is a length 5 vector containing the parameters of the modified Givens rotation:
/// param(0) = flag, param(1) = h11, param(2) = h21, param(3) = h12, param(4) = h22
///
/// \tparam MemberType: TeamPolicy member type
template <typename MemberType>
struct TeamRotm {
  /// \brief Invokes the TeamRotm functor. A nested parallel_for with TeamThreadRange is used.
  /// \tparam XViewType: Input/output type for the vector x, needs to be a 1D view
  /// \tparam YViewType: Input/output type for the vector y, needs to be a 1D view
  /// \tparam ParamViewType: Input type for the param vector, needs to be a 1D view of length 5
  ///
  /// \param[in,out] x: x is a length n vector, a rank 1 view
  /// \param[in,out] y: y is a length n vector, a rank 1 view
  /// \param[in] param: param is a length 5 vector, a rank 1 view, containing the parameters of the modified Givens
  /// rotation
  template <typename XViewType, typename YViewType, typename ParamViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &x, const YViewType &y,
                                           const ParamViewType &param);
};

/// \brief TeamVector Batched Rotm:
/// Applies the modified Givens rotation to vectors x and y:
///  x(i) := h11*x(i) + h12*y(i)
///  y(i) := h21*x(i) + h22*y(i)
///
/// The matrix H is given by
/// flag == -1.0  flag ==  0.0  flag ==  1.0  flag == -2.0
/// [[h11, h12],  [[1, h12],    [[h11, 1],    [[1, 0],
///  [h21, h22]]   [h21, 1]]     [[-1, h12]]   [[0, 1]]
///
/// param is a length 5 vector containing the parameters of the modified Givens rotation:
/// param(0) = flag, param(1) = h11, param(2) = h21, param(3) = h12, param(4) = h22
///
/// \tparam MemberType: TeamPolicy member type
template <typename MemberType>
struct TeamVectorRotm {
  /// \brief Invokes the TeamVectorRotm functor. A nested parallel_for with TeamVectorRange is used.
  /// \tparam MemberType: TeamPolicy member type
  /// \tparam XViewType: Input/output type for the vector x, needs to be a 1D view
  /// \tparam YViewType: Input/output type for the vector y, needs to be a 1D view
  /// \tparam ParamViewType: Input type for the param vector, needs to be a 1D view of length 5
  ///
  /// \param[in,out] x: x is a length n vector, a rank 1 view
  /// \param[in,out] y: y is a length n vector, a rank 1 view
  /// \param[in] param: param is a length 5 vector, a rank 1 view, containing the parameters of the modified Givens
  /// rotation
  template <typename XViewType, typename YViewType, typename ParamViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &x, const YViewType &y,
                                           const ParamViewType &param);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Rotm_Impl.hpp"

#endif  // KOKKOSBATCHED_ROTM_HPP_
