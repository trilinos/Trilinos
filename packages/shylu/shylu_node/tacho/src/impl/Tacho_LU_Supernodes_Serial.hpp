// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_LU_SUPERNODES_SERIAL_HPP__
#define __TACHO_LU_SUPERNODES_SERIAL_HPP__

/// \file Tacho_LU_Supernodes.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

#include "Tacho_Lapack_External.hpp"
#include "Tacho_Lapack_Team.hpp"

#include "Tacho_Blas_External.hpp"
#include "Tacho_Blas_Team.hpp"

#include "Tacho_ApplyPivots.hpp"
#include "Tacho_ApplyPivots_Internal.hpp"

#include "Tacho_LU.hpp"
#include "Tacho_LU_External.hpp"
#include "Tacho_LU_Internal.hpp"

#include "Tacho_Trsm.hpp"
#include "Tacho_Trsm_External.hpp"
#include "Tacho_Trsm_Internal.hpp"

#include "Tacho_Gemm.hpp"
#include "Tacho_Gemm_External.hpp"
#include "Tacho_Gemm_Internal.hpp"

#include "Tacho_Trsv.hpp"
#include "Tacho_Trsv_External.hpp"
#include "Tacho_Trsv_Internal.hpp"

#include "Tacho_Gemv.hpp"
#include "Tacho_Gemv_External.hpp"
#include "Tacho_Gemv_Internal.hpp"

namespace Tacho {

template <> struct LU_Supernodes<Algo::Workflow::Serial> {
  template <typename MemberType, typename SupernodeInfoType>
  KOKKOS_INLINE_FUNCTION static int
  factorize(MemberType &member, const SupernodeInfoType &info, const typename SupernodeInfoType::ordinal_type_array &P,
            const typename SupernodeInfoType::value_type_matrix &ABR, const ordinal_type sid) {
    using supernode_info_type = SupernodeInfoType;
    using value_type = typename supernode_info_type::value_type;
    using value_type_matrix = typename supernode_info_type::value_type_matrix;
    using range_type = typename supernode_info_type::range_type;

    using LU_AlgoType = typename LU_Algorithm::type;
    using TrsmAlgoType = typename TrsmAlgorithm::type;
    using GemmAlgoType = typename GemmAlgorithm::type;

    // get current supernode
    const auto &s = info.supernodes(sid);

    // panel (s.m x s.n) is divided into ATL (m x m) and ATR (m x n)
    const ordinal_type m = s.m, n = s.n, n_m = s.n - s.m;

    // m and n are available, then factorize the supernode block
    if (m > 0) {
      /// LU factorize ATL
      UnmanagedViewType<value_type_matrix> AT(s.u_buf, m, n);

      LU<LU_AlgoType>::invoke(member, AT, P);
      LU<LU_AlgoType>::modify(member, m, P);

      if (n > 0) {
        const value_type one(1), minus_one(-1), zero(0);

        UnmanagedViewType<value_type_matrix> ATL(s.u_buf, m, m);
        UnmanagedViewType<value_type_matrix> AL(s.l_buf, n, m);
        const auto ABL = Kokkos::subview(AL, range_type(m, n), Kokkos::ALL());

        Trsm<Side::Right, Uplo::Upper, Trans::NoTranspose, TrsmAlgoType>::invoke(member, Diag::NonUnit(), one, ATL,
                                                                                 ABL);

        TACHO_TEST_FOR_ABORT(static_cast<ordinal_type>(ABR.extent(0)) != n_m ||
                                 static_cast<ordinal_type>(ABR.extent(1)) != n_m,
                             "ABR dimension does not match to supernodes");

        UnmanagedViewType<value_type_matrix> ATR(s.u_buf + ATL.span(), m, n_m);
        Gemm<Trans::NoTranspose, Trans::NoTranspose, GemmAlgoType>::invoke(member, minus_one, ABL, ATR, zero, ABR);
      }
    }
    return 0;
  }

  template <typename MemberType, typename SupernodeInfoType>
  KOKKOS_INLINE_FUNCTION static int solve_lower(MemberType &member, const SupernodeInfoType &info,
                                                const typename SupernodeInfoType::ordinal_type_array &P,
                                                const typename SupernodeInfoType::value_type_matrix &xB,
                                                const ordinal_type sid) {
    using supernode_info_type = SupernodeInfoType;
    using value_type = typename supernode_info_type::value_type;
    using value_type_matrix = typename supernode_info_type::value_type_matrix;
    using ordinal_type_array = typename supernode_info_type::ordinal_type_array;
    using range_type = typename supernode_info_type::range_type;

    const auto &s = info.supernodes(sid);

    using TrsvAlgoType = typename TrsvAlgorithm::type;
    using GemvAlgoType = typename GemvAlgorithm::type;

    // panel is divided into diagonal and interface block
    const ordinal_type m = s.m, n = s.n, n_m = s.n - s.m; //, nrhs = info.x.extent(1);

    // m and n are available, then factorize the supernode block
    if (m > 0) {
      const value_type minus_one(-1), zero(0);
      const ordinal_type offm = s.row_begin;
      UnmanagedViewType<value_type_matrix> ATL(s.u_buf, m, m);
      const auto xT = Kokkos::subview(info.x, range_type(offm, offm + m), Kokkos::ALL());
      const ConstUnmanagedViewType<ordinal_type_array> fpiv(P.data() + m, m);

      ApplyPivots<PivotMode::Flame, Side::Left, Direct::Forward, Algo::Internal> /// row inter-change
          ::invoke(member, fpiv, xT);

      Trsv<Uplo::Lower, Trans::NoTranspose, TrsvAlgoType>::invoke(member, Diag::Unit(), ATL, xT);

      if (n_m > 0) {
        UnmanagedViewType<value_type_matrix> AL(s.l_buf, n, m);
        const auto ABL = Kokkos::subview(AL, range_type(m, n), Kokkos::ALL());
        Gemv<Trans::NoTranspose, GemvAlgoType>::invoke(member, minus_one, ABL, xT, zero, xB);
      }
    }
    return 0;
  }

  template <typename MemberType, typename SupernodeInfoType>
  KOKKOS_INLINE_FUNCTION static int solve_upper(MemberType &member, const SupernodeInfoType &info,
                                                const typename SupernodeInfoType::value_type_matrix &xB,
                                                const ordinal_type sid) {
    using supernode_info_type = SupernodeInfoType;

    using value_type = typename supernode_info_type::value_type;
    using value_type_matrix = typename supernode_info_type::value_type_matrix;
    using range_type = typename supernode_info_type::range_type;

    using TrsvAlgoType = typename TrsvAlgorithm::type;
    using GemvAlgoType = typename GemvAlgorithm::type;

    // get current supernode
    const auto &s = info.supernodes(sid);

    // panel is divided into diagonal and interface block
    const ordinal_type m = s.m, n_m = s.n - s.m;

    // m and n are available, then factorize the supernode block
    if (m > 0) {
      const value_type minus_one(-1), one(1);
      const UnmanagedViewType<value_type_matrix> ATL(s.u_buf, m, m);

      const ordinal_type offm = s.row_begin;
      const auto xT = Kokkos::subview(info.x, range_type(offm, offm + m), Kokkos::ALL());

      if (n_m > 0) {
        const UnmanagedViewType<value_type_matrix> ATR(s.u_buf + ATL.span(), m, n_m); // ptr += m*n;
        Gemv<Trans::NoTranspose, GemvAlgoType>::invoke(member, minus_one, ATR, xB, one, xT);
      }
      Trsv<Uplo::Upper, Trans::NoTranspose, TrsvAlgoType>::invoke(member, Diag::NonUnit(), ATL, xT);
    }
    return 0;
  }

  template <typename MemberType, typename SupernodeInfoType>
  KOKKOS_INLINE_FUNCTION static int
  factorize_recursive_serial(MemberType &member, const SupernodeInfoType &info, const ordinal_type sid,
                             const bool final, typename SupernodeInfoType::ordinal_type_array::pointer_type piv,
                             typename SupernodeInfoType::value_type_array::pointer_type buf, const size_type bufsize) {
    using supernode_info_type = SupernodeInfoType;
    using value_type = typename supernode_info_type::value_type;
    using value_type_matrix = typename supernode_info_type::value_type_matrix;
    using ordinal_type_array = typename supernode_info_type::ordinal_type_array;

    const auto &s = info.supernodes(sid);
    if (final) {
      // serial recursion
      for (ordinal_type i = 0; i < s.nchildren; ++i)
        factorize_recursive_serial(member, info, s.children[i], final, piv, buf, bufsize);
    }

    {
      const ordinal_type m = s.m;
      const ordinal_type rbeg = s.row_begin;
      UnmanagedViewType<ordinal_type_array> ipiv(piv + rbeg * 4, 4 * m);

      const ordinal_type n = s.n - s.m;
      const size_type bufsize_required = n * n * sizeof(value_type);
      TACHO_TEST_FOR_ABORT(bufsize < bufsize_required, "bufsize is smaller than required");
      value_type *bufptr = buf;
      UnmanagedViewType<value_type_matrix> ABR(bufptr, n, n);

      LU_Supernodes<Algo::Workflow::Serial>::factorize(member, info, ipiv, ABR, sid);

      constexpr bool update_lower = true;
      CholSupernodes<Algo::Workflow::Serial>::update(member, info, ABR, sid, bufsize - ABR.span() * sizeof(value_type),
                                                     (void *)((value_type *)buf + ABR.span()), update_lower);
    }
    return 0;
  }

  template <typename MemberType, typename SupernodeInfoType>
  KOKKOS_INLINE_FUNCTION static int
  solve_lower_recursive_serial(MemberType &member, const SupernodeInfoType &info, const ordinal_type sid,
                               const bool final, typename SupernodeInfoType::ordinal_type_array::pointer_type piv,
                               typename SupernodeInfoType::value_type_array::pointer_type buf,
                               const size_type bufsize) {
    using supernode_info_type = SupernodeInfoType;

    using value_type = typename supernode_info_type::value_type;
    using value_type_matrix = typename supernode_info_type::value_type_matrix;
    using ordinal_type_array = typename supernode_info_type::ordinal_type_array;

    const auto &s = info.supernodes(sid);

    if (final) {
      // serial recursion
      for (ordinal_type i = 0; i < s.nchildren; ++i)
        solve_lower_recursive_serial(member, info, s.children[i], final, piv, buf, bufsize);
    }

    {
      const ordinal_type m = s.m;
      const ordinal_type rbeg = s.row_begin;
      UnmanagedViewType<ordinal_type_array> P(piv + rbeg * 4, 4 * m);

      const ordinal_type n = s.n - s.m;
      const ordinal_type nrhs = info.x.extent(1);
      const size_type bufsize_required = n * nrhs * sizeof(value_type);

      TACHO_TEST_FOR_ABORT(bufsize < bufsize_required, "bufsize is smaller than required");

      UnmanagedViewType<value_type_matrix> xB((value_type *)buf, n, nrhs);

      LU_Supernodes<Algo::Workflow::Serial>::solve_lower(member, info, P, xB, sid);

      CholSupernodes<Algo::Workflow::Serial>::update_solve_lower(member, info, xB, sid);
    }
    return 0;
  }

  template <typename MemberType, typename SupernodeInfoType>
  KOKKOS_INLINE_FUNCTION static int
  solve_upper_recursive_serial(MemberType &member, const SupernodeInfoType &info, const ordinal_type sid,
                               const bool final, typename SupernodeInfoType::ordinal_type_array::pointer_type piv,
                               typename SupernodeInfoType::value_type_array::pointer_type buf,
                               const ordinal_type bufsize) {
    using supernode_info_type = SupernodeInfoType;
    using value_type = typename supernode_info_type::value_type;
    using value_type_matrix = typename supernode_info_type::value_type_matrix;

    const auto &s = info.supernodes(sid);
    {
      const ordinal_type n = s.n - s.m;
      const ordinal_type nrhs = info.x.extent(1);
      const ordinal_type bufsize_required = n * nrhs * sizeof(value_type);

      TACHO_TEST_FOR_ABORT(bufsize < bufsize_required, "bufsize is smaller than required");

      UnmanagedViewType<value_type_matrix> xB((value_type *)buf, n, nrhs);

      CholSupernodes<Algo::Workflow::Serial>::update_solve_upper(member, info, xB, sid);

      LU_Supernodes<Algo::Workflow::Serial>::solve_upper(member, info, xB, sid);
    }

    if (final) {
      // serial recursion
      for (ordinal_type i = 0; i < s.nchildren; ++i)
        solve_upper_recursive_serial(member, info, s.children[i], final, piv, buf, bufsize);
    }
    return 0;
  }
};
} // namespace Tacho

#endif
