// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef __KOKKOSBATCHED_SOLVE_UTV_TEAMVECTOR_IMPL_COMPADRE_HPP__
#define __KOKKOSBATCHED_SOLVE_UTV_TEAMVECTOR_IMPL_COMPADRE_HPP__

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_SolveUTV_TeamVector_Internal_Compadre.hpp"

namespace KokkosBatched {

    ///
    /// TeamVector Impl
    /// ===============
    template<typename MemberType>
    struct TeamVectorSolveUTVCompadre<MemberType,Algo::UTV::Unblocked> {
        template<typename UViewType,
             typename TViewType,
             typename VViewType,
             typename pViewType,
             typename BViewType,
             typename XViewType,
             typename wViewType>
        KOKKOS_INLINE_FUNCTION
        static int
        invoke(const MemberType &member, 
             const int matrix_rank,
             const int m,
             const int n,
             const int nrhs,
             const UViewType &U,
             const TViewType &T,
             const VViewType &V,
             const pViewType &p,
             const BViewType &B,
             const XViewType &X,
             const wViewType &w_a,
             const wViewType &w_b,
             const bool implicit_RHS) {
              TeamVectorSolveUTV_Internal_Compadre::
                invoke(member,
                   matrix_rank, m, n, nrhs,
                   U.data(), U.stride(0), U.stride(1),
                   T.data(), T.stride(0), T.stride(1),
                   V.data(), V.stride(0), V.stride(1),
                   p.data(), p.stride(0),
                   B.data(), B.stride(0), B.stride(1),
                   X.data(), X.stride(0), X.stride(1),
                   w_a.data(), w_b.data(),
                   implicit_RHS);
                return 0;
        }
    };
}



#endif
