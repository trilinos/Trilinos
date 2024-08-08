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
#ifndef __KOKKOSBATCHED_UTV_DECL_HPP__
#define __KOKKOSBATCHED_UTV_DECL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// For given A, it performs UTV factorization i.e., UTV = A P^T
/// - input:
///   - A is m x m real matrix
///   - p is m integer vector
///   - U is m x m real matrix
///   - V is m x m real matrix
///   - w is 3*m real vector workspace (contiguous)
/// - output:
///   - A is overwritten as lower triangular matrix_rank x matrix_rank real
///   matrix
///   - P^T includes pivot indicies (note that this is different from
///   permutation indicies)
///   - U is left orthogonal matrix m x matrix_rank
///   - V is right orthogonal matrix matrix_rank x m
///
/// When A is a full rank i.e., matrix_rank == m, this only compute a QR with
/// column pivoting
/// - output:
///   - A is overwritten as upper triangular matrix
///   - P^T includes pivot indicies (note that this is different from
///   permutation indicies)
///   - U is an orthogonal matrix m x m
///   - V is not touched
///
/// For the solution of a rank-deficient problem, it is recommended to use
/// SolveUTV.
///

///
/// TeamVector UTV
///

template <typename MemberType, typename ArgAlgo>
struct TeamVectorUTV {
  template <typename AViewType, typename pViewType, typename UViewType, typename VViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const pViewType &p,
                                           const UViewType &U, const VViewType &V, const wViewType &w,
                                           int &matrix_rank);
};

}  // namespace KokkosBatched

#include "KokkosBatched_UTV_TeamVector_Impl.hpp"

#endif
