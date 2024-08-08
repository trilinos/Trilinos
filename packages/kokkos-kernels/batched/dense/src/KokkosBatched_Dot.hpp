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
#ifndef __KOKKOSBATCHED_DOT_HPP__
#define __KOKKOSBATCHED_DOT_HPP__

/// \author Kim Liegeois (knliege@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

/// \brief Serial Batched DOT:
///
/// Depending on the ArgTrans template, the dot product is
/// row-based (ArgTrans == Trans::NoTranspose):
///
///   dot_l <- (x_l:, y_l:) for all l = 1, ..., N
/// where:
///   * N is the second dimension of X.
///
/// Or column-based:
///   dot_l <- (x_:l, y_:l) for all l = 1, ..., n
/// where:
///   * n is the second dimension of X.
///
/// \tparam ArgTrans: type of dot product (Trans::NoTranspose by default)
/// \tparam XViewType: Input type for X, needs to be a 2D view
/// \tparam YViewType: Input type for Y, needs to be a 2D view
/// \tparam alphaViewType: Input type for alpha, needs to be a 1D view
///
/// \param X [in]: Input vector X, a rank 2 view
/// \param Y [in]: Input vector Y, a rank 2 view
/// \param dot [out]: Computed dot product, a rank 1 view
///
/// No nested parallel_for is used inside of the function.
///

template <typename ArgTrans = Trans::NoTranspose>
struct SerialDot {
  template <typename XViewType, typename YViewType, typename NormViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const XViewType &X, const YViewType &Y, const NormViewType &dot);
};

/// \brief Team Batched DOT:
///
/// Depending on the ArgTrans template, the dot product is
/// row-based (ArgTrans == Trans::NoTranspose):
///
///   dot_l <- (x_l:, y_l:) for all l = 1, ..., N
/// where:
///   * N is the second dimension of X.
///
/// Or column-based:
///   dot_l <- (x_:l, y_:l) for all l = 1, ..., n
/// where:
///   * n is the second dimension of X.
///
/// \tparam ArgTrans: type of dot product (Trans::NoTranspose by default)
/// \tparam XViewType: Input type for X, needs to be a 2D view
/// \tparam YViewType: Input type for Y, needs to be a 2D view
/// \tparam alphaViewType: Input type for alpha, needs to be a 1D view
///
/// \param X [in]: Input vector X, a rank 2 view
/// \param Y [in]: Input vector Y, a rank 2 view
/// \param dot [out]: Computed dot product, a rank 1 view
///
/// A nested parallel_for with TeamThreadRange is used.
///

template <typename MemberType, typename ArgTrans = Trans::NoTranspose>
struct TeamDot {
  template <typename XViewType, typename YViewType, typename NormViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &X, const YViewType &Y,
                                           const NormViewType &dot);
};

/// \brief TeamVector Batched DOT:
///
/// Depending on the ArgTrans template, the dot product is
/// row-based (ArgTrans == Trans::NoTranspose):
///
///   dot_l <- (x_l:, y_l:) for all l = 1, ..., N
/// where:
///   * N is the second dimension of X.
///
/// Or column-based:
///   dot_l <- (x_:l, y_:l) for all l = 1, ..., n
/// where:
///   * n is the second dimension of X.
///
/// \tparam ArgTrans: type of dot product (Trans::NoTranspose by default)
/// \tparam XViewType: Input type for X, needs to be a 2D view
/// \tparam YViewType: Input type for Y, needs to be a 2D view
/// \tparam alphaViewType: Input type for alpha, needs to be a 1D view
///
/// \param X [in]: Input vector X, a rank 2 view
/// \param Y [in]: Input vector Y, a rank 2 view
/// \param dot [out]: Computed dot product, a rank 1 view
///
/// Two nested parallel_for with both TeamThreadRange and ThreadVectorRange
/// (or one with TeamVectorRange) are used inside.
///

template <typename MemberType, typename ArgTrans = Trans::NoTranspose>
struct TeamVectorDot {
  template <typename XViewType, typename YViewType, typename NormViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &X, const YViewType &Y,
                                           const NormViewType &dot);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Dot_Internal.hpp"

#endif
