//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#ifndef __KOKKOSBATCHED_AXPY_HPP__
#define __KOKKOSBATCHED_AXPY_HPP__

/// \author Kim Liegeois (knliege@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

/// \brief Serial Batched AXPY:
///   y_l <- alpha_l * x_l + y_l for all l = 1, ..., N
/// where:
///   * N is the number of vectors,
///   * x_1, ..., x_N are the N input vectors,
///   * y_1, ..., y_N are the N output vectors,
///   * alpha_1, ..., alpha_N are N scaling factors for x_1, ..., x_N.
///
/// \tparam XViewType: Input type for X, needs to be a 2D view
/// \tparam YViewType: Input type for Y, needs to be a 2D view
/// \tparam alphaViewType: Input type for alpha, needs to be a 1D view
///
/// \param alpha [in]: input coefficient for X, a rank 1 view
/// \param X [in]: Input vector X, a rank 2 view
/// \param Y [in/out]: Output vector Y, a rank 2 view
///
/// No nested parallel_for is used inside of the function.
///

struct SerialAxpy {
  template <typename XViewType, typename YViewType, typename alphaViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const alphaViewType &alpha, const XViewType &X, const YViewType &Y);
};

/// \brief Team Batched AXPY:
///   y_l <- alpha_l * x_l + y_l for all l = 1, ..., N
/// where:
///   * N is the number of vectors,
///   * x_1, ..., x_N are the N input vectors,
///   * y_1, ..., y_N are the N output vectors,
///   * alpha_1, ..., alpha_N are N scaling factors for x_1, ..., x_N.
///
/// \tparam XViewType: Input type for X, needs to be a 2D view
/// \tparam YViewType: Input type for Y, needs to be a 2D view
/// \tparam alphaViewType: Input type for alpha, needs to be a 1D view
///
/// \param member [in]: TeamPolicy member
/// \param alpha [in]: input coefficient for X, a rank 1 view
/// \param X [in]: Input vector X, a rank 2 view
/// \param Y [in/out]: Output vector Y, a rank 2 view
///
/// A nested parallel_for with TeamThreadRange is used.
///

template <typename MemberType>
struct TeamAxpy {
  template <typename XViewType, typename YViewType, typename alphaViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const alphaViewType &alpha, const XViewType &X,
                                           const YViewType &Y);
};

/// \brief TeamVector Batched AXPY:
///   y_l <- alpha_l * x_l + y_l for all l = 1, ..., N
/// where:
///   * N is the number of vectors,
///   * x_1, ..., x_N are the N input vectors,
///   * y_1, ..., y_N are the N output vectors,
///   * alpha_1, ..., alpha_N are N scaling factors for x_1, ..., x_N.
///
/// \tparam XViewType: Input type for X, needs to be a 2D view
/// \tparam YViewType: Input type for Y, needs to be a 2D view
/// \tparam alphaViewType: Input type for alpha, needs to be a 1D view
///
/// \param member [in]: TeamPolicy member
/// \param alpha [in]: input coefficient for X, a rank 1 view
/// \param X [in]: Input vector X, a rank 2 view
/// \param Y [in/out]: Output vector Y, a rank 2 view
///
/// Two nested parallel_for with both TeamThreadRange and ThreadVectorRange
/// (or one with TeamVectorRange) are used inside.
///

template <typename MemberType>
struct TeamVectorAxpy {
  template <typename XViewType, typename YViewType, typename alphaViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const alphaViewType &alpha, const XViewType &X,
                                           const YViewType &Y);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Axpy_Impl.hpp"

#endif
