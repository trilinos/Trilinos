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
#ifndef __KOKKOSBATCHED_CG_HPP__
#define __KOKKOSBATCHED_CG_HPP__

/// \author Kim Liegeois (knliege@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

/// \brief Batched CG: Selective Interface
///
/// \tparam OperatorType: The type of the operator of the system
/// \tparam VectorViewType: Input type for the right-hand side and the solution,
/// needs to be a 2D view
///
/// \param member [in]: TeamPolicy member
/// \param A [in]: batched operator (can be a batched matrix or a (left or right
/// or both) preconditioned batched matrix) \param B [in]: right-hand side, a
/// rank 2 view \param X [in/out]: initial guess and solution, a rank 2 view
/// \param handle [in]: a handle which provides different information such as
/// the tolerance or the maximal number of iterations of the solver.

#include "KokkosBatched_Krylov_Handle.hpp"
#include "KokkosBatched_CG_Team_Impl.hpp"
#include "KokkosBatched_CG_TeamVector_Impl.hpp"

namespace KokkosBatched {

template <typename MemberType, typename ArgMode>
struct CG {
  template <typename OperatorType, typename VectorViewType, typename KrylovHandleType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const OperatorType &A, const VectorViewType &B,
                                           const VectorViewType &X, const KrylovHandleType &handle) {
    int status = 0;
    if (std::is_same<ArgMode, Mode::Team>::value) {
      status = TeamCG<MemberType>::template invoke<OperatorType, VectorViewType>(member, A, B, X, handle);
    } else if (std::is_same<ArgMode, Mode::TeamVector>::value) {
      status = TeamVectorCG<MemberType>::template invoke<OperatorType, VectorViewType>(member, A, B, X, handle);
    }
    return status;
  }
};

}  // namespace KokkosBatched
#endif
