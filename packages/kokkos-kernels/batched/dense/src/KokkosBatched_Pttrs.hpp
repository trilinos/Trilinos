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
#ifndef KOKKOSBATCHED_PTTRS_HPP_
#define KOKKOSBATCHED_PTTRS_HPP_

#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Serial Batched Pttrs:
/// Solve Ab_l x_l = b_l for all l = 0, ..., N
///   using the factorization A = U**H * D * U or A = L * D * L**H computed by
///   Pttrf.
///
/// \tparam DViewType: Input type for the a diagonal matrix, needs to be a 1D
/// view
/// \tparam EViewType: Input type for the a upper/lower diagonal matrix,
/// needs to be a 1D view
/// \tparam BViewType: Input type for the right-hand side and the solution,
/// needs to be a 1D view
///
/// \param d [in]: n diagonal elements of the diagonal matrix D
/// \param e [in]: n-1 upper/lower diagonal elements of the diagonal matrix E
/// \param b [inout]: right-hand side and the solution, a rank 1 view
///
/// No nested parallel_for is used inside of the function.
///

template <typename ArgUplo, typename ArgAlgo>
struct SerialPttrs {
  template <typename DViewType, typename EViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const DViewType &d, const EViewType &e, const BViewType &b);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Pttrs_Serial_Impl.hpp"

#endif  // KOKKOSBATCHED_PTTRS_HPP_
