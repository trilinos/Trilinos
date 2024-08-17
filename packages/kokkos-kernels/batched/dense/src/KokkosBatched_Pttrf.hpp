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
#ifndef KOKKOSBATCHED_PTTRF_HPP_
#define KOKKOSBATCHED_PTTRF_HPP_

#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Serial Batched Pttrf:
/// Compute the Cholesky factorization L*D*L**T (or L*D*L**H) of a real
/// symmetric (or complex Hermitian) positive definite tridiagonal matrix A_l
/// for all l = 0, ..., N
///
/// \tparam DViewType: Input type for the a diagonal matrix, needs to be a 1D
/// view
/// \tparam EViewType: Input type for the a upper/lower diagonal matrix,
/// needs to be a 1D view
///
/// \param d [inout]: n diagonal elements of the diagonal matrix D
/// \param e [inout]: n-1 upper/lower diagonal elements of the diagonal matrix E
///
/// No nested parallel_for is used inside of the function.
///

template <typename ArgAlgo>
struct SerialPttrf {
  template <typename DViewType, typename EViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const DViewType &d, const EViewType &e);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Pttrf_Serial_Impl.hpp"

#endif  // KOKKOSBATCHED_PTTRF_HPP_
