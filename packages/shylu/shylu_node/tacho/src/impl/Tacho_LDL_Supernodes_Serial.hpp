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
#ifndef __TACHO_LDL_SUPERNODES_SERIAL_HPP__
#define __TACHO_LDL_SUPERNODES_SERIAL_HPP__

/// \file Tacho_LDL_Supernodes.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_CholSupernodes_Serial.hpp"

#include "Tacho_Symmetrize.hpp"
#include "Tacho_Symmetrize_Internal.hpp"

#include "Tacho_ApplyPivots.hpp"
#include "Tacho_ApplyPivots_Internal.hpp"

#include "Tacho_Copy.hpp"
#include "Tacho_Copy_Internal.hpp"

#include "Tacho_Scale2x2_BlockInverseDiagonals.hpp"
#include "Tacho_Scale2x2_BlockInverseDiagonals_Internal.hpp"

#include "Tacho_LDL.hpp"
#include "Tacho_LDL_External.hpp"
#include "Tacho_LDL_Internal.hpp"

#include "Tacho_GemmTriangular.hpp"
#include "Tacho_GemmTriangular_External.hpp"
#include "Tacho_GemmTriangular_Internal.hpp"

namespace Tacho {

template <> struct LDL_Supernodes<Algo::Workflow::Serial> {
  template <typename MemberType, typename SupernodeInfoType>
  KOKKOS_INLINE_FUNCTION static int
  factorize(MemberType &member, const SupernodeInfoType &info, const typename SupernodeInfoType::ordinal_type_array &P,
            const typename SupernodeInfoType::value_type_matrix &D,
            const typename SupernodeInfoType::value_type_array &W,
            const typename SupernodeInfoType::value_type_matrix &ABR, const ordinal_type sid) {
    using supernode_info_type = SupernodeInfoType;
    using value_type = typename supernode_info_type::value_type;
    using value_type_matrix = typename supernode_info_type::value_type_matrix;
    using ordinal_type_array = typename supernode_info_type::ordinal_type_array;

    // algorithm choice
    using LDL_AlgoType = typename LDL_Algorithm::type;
    using TrsmAlgoType = typename TrsmAlgorithm::type;
    using GemmAlgoType = typename GemmAlgorithm::type;

    // get current supernode
    const auto &s = info.supernodes(sid);

    // get panel pointer
    value_type *ptr = s.u_buf;

    // panel (s.m x s.n) is divided into ATL (m x m) and ATR (m x n)
    const ordinal_type m = s.m, n = s.n - s.m;

    // m and n are available, then factorize the supernode block
    if (m > 0) {
      /// LDL factorize ATL, extract diag, symmetrize ATL with unit diagonals
      UnmanagedViewType<value_type_matrix> ATL(ptr, m, m);
      ptr += m * m;

      Symmetrize<Uplo::Upper, Algo::Internal>::invoke(member, ATL);

      LDL<Uplo::Lower, LDL_AlgoType>::invoke(member, ATL, P, W);
      LDL<Uplo::Lower, LDL_AlgoType>::modify(member, ATL, P, D);

      if (n > 0) {
        const value_type one(1), zero(0);
        UnmanagedViewType<value_type_matrix> ATR(ptr, m, n);
        ptr += m * n;
        UnmanagedViewType<value_type_matrix> STR(W.data(), m, n);

        auto fpiv = ordinal_type_array(P.data() + m, m);
        ApplyPivots<PivotMode::Flame, Side::Left, Direct::Forward, Algo::Internal> /// row inter-change
            ::invoke(member, fpiv, ATR);
        Trsm<Side::Left, Uplo::Lower, Trans::NoTranspose, TrsmAlgoType>::invoke(member, Diag::Unit(), one, ATL, ATR);
        Copy<Algo::Internal>::invoke(member, STR, ATR);
        Scale2x2_BlockInverseDiagonals<Side::Left, Algo::Internal> /// row scaling
            ::invoke(member, P, D, ATR);

        TACHO_TEST_FOR_ABORT(static_cast<ordinal_type>(ABR.extent(0)) != n ||
                                 static_cast<ordinal_type>(ABR.extent(1)) != n,
                             "ABR dimension does not match to supernodes");
        GemmTriangular<Trans::Transpose, Trans::NoTranspose, Uplo::Upper, GemmAlgoType>::invoke(member, -one, ATR, STR,
                                                                                                zero, ABR);
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

    using range_type = Kokkos::pair<ordinal_type, ordinal_type>;

    const auto &s = info.supernodes(sid);

    using TrsvAlgoType = typename TrsvAlgorithm::type;
    using GemvAlgoType = typename GemvAlgorithm::type;

    // get panel pointer
    value_type *ptr = s.u_buf;

    // panel is divided into diagonal and interface block
    const ordinal_type m = s.m, n = s.n - s.m; //, nrhs = info.x.extent(1);

    // m and n are available, then factorize the supernode block
    if (m > 0) {
      const value_type one(1), zero(0);
      const ordinal_type offm = s.row_begin;
      UnmanagedViewType<value_type_matrix> AL(ptr, m, m);
      ptr += m * m;
      const auto xT = Kokkos::subview(info.x, range_type(offm, offm + m), Kokkos::ALL());
      const auto fpiv = ordinal_type_array(P.data() + m, m);

      ApplyPivots<PivotMode::Flame, Side::Left, Direct::Forward, Algo::Internal> /// row inter-change
          ::invoke(member, fpiv, xT);

      Trsv<Uplo::Lower, Trans::NoTranspose, TrsvAlgoType>::invoke(member, Diag::Unit(), AL, xT);

      if (n > 0) {
        UnmanagedViewType<value_type_matrix> AR(ptr, m, n); // ptr += m*n;
        Gemv<Trans::Transpose, GemvAlgoType>::invoke(member, -one, AR, xT, zero, xB);
      }
    }
    return 0;
  }

  template <typename MemberType, typename SupernodeInfoType>
  KOKKOS_INLINE_FUNCTION static int solve_upper(MemberType &member, const SupernodeInfoType &info,
                                                const typename SupernodeInfoType::ordinal_type_array &P,
                                                const typename SupernodeInfoType::value_type_matrix &D,
                                                const typename SupernodeInfoType::value_type_matrix &xB,
                                                const ordinal_type sid) {
    using supernode_info_type = SupernodeInfoType;

    using value_type = typename supernode_info_type::value_type;
    using value_type_matrix = typename supernode_info_type::value_type_matrix;
    using ordinal_type_array = typename supernode_info_type::ordinal_type_array;

    using range_type = Kokkos::pair<ordinal_type, ordinal_type>;

    using GemvAlgoType = typename GemvAlgorithm::type;
    using TrsvAlgoType = typename TrsvAlgorithm::type;

    // get current supernode
    const auto &s = info.supernodes(sid);

    // get supernode panel pointer
    value_type *ptr = s.u_buf;

    // panel is divided into diagonal and interface block
    const ordinal_type m = s.m, n = s.n - s.m; //, nrhs = info.x.extent(1);

    // m and n are available, then factorize the supernode block
    if (m > 0) {
      const value_type one(1);
      const UnmanagedViewType<value_type_matrix> AL(ptr, m, m);
      ptr += m * m;

      const ordinal_type offm = s.row_begin;
      const auto xT = Kokkos::subview(info.x, range_type(offm, offm + m), Kokkos::ALL());
      const auto fpiv = ordinal_type_array(P.data() + m, m);

      Scale2x2_BlockInverseDiagonals<Side::Left, Algo::Internal> /// row scaling
          ::invoke(member, P, D, xT);

      if (n > 0) {
        const UnmanagedViewType<value_type_matrix> AR(ptr, m, n); // ptr += m*n;
        Gemv<Trans::NoTranspose, GemvAlgoType>::invoke(member, -one, AR, xB, one, xT);
      }
      Trsv<Uplo::Lower, Trans::Transpose, TrsvAlgoType>::invoke(member, Diag::Unit(), AL, xT);
      ApplyPivots<PivotMode::Flame, Side::Left, Direct::Backward, Algo::Internal> /// row inter-change
          ::invoke(member, fpiv, xT);
    }
    return 0;
  }

  template <typename MemberType, typename SupernodeInfoType>
  KOKKOS_INLINE_FUNCTION static int
  factorize_recursive_serial(MemberType &member, const SupernodeInfoType &info, const ordinal_type sid,
                             const bool final, typename SupernodeInfoType::ordinal_type_array::pointer_type piv,
                             typename SupernodeInfoType::value_type_array::pointer_type diag,
                             typename SupernodeInfoType::value_type_array::pointer_type buf, const size_type bufsize) {
    using supernode_info_type = SupernodeInfoType;
    using value_type = typename supernode_info_type::value_type;
    using value_type_array = typename supernode_info_type::value_type_array;
    using value_type_matrix = typename supernode_info_type::value_type_matrix;
    using ordinal_type_array = typename supernode_info_type::ordinal_type_array;

    const auto &s = info.supernodes(sid);
    if (final) {
      // serial recursion
      for (ordinal_type i = 0; i < s.nchildren; ++i)
        factorize_recursive_serial(member, info, s.children[i], final, piv, diag, buf, bufsize);
    }

    {
      const ordinal_type m = s.m;
      const ordinal_type rbeg = s.row_begin;
      UnmanagedViewType<ordinal_type_array> ipiv(piv + rbeg * 4, 4 * m);
      UnmanagedViewType<value_type_matrix> dblk(diag + rbeg * 2, m, 2);

      const ordinal_type n = s.n - s.m;

      const ordinal_type mm = m < 32 ? m : 32;
      const ordinal_type mn = mm > n ? mm : n;

      const size_type bufsize_required = (n * n + m * mn) * sizeof(value_type);
      TACHO_TEST_FOR_ABORT(bufsize < bufsize_required, "bufsize is smaller than required");
      value_type *bufptr = buf;
      UnmanagedViewType<value_type_matrix> ABR(bufptr, n, n);
      bufptr += ABR.span();
      UnmanagedViewType<value_type_array> w(bufptr, m * mn);
      bufptr += w.span();

      LDL_Supernodes<Algo::Workflow::Serial>::factorize(member, info, ipiv, dblk, w, ABR, sid);

      /// assembly is same
      CholSupernodes<Algo::Workflow::Serial>::update(member, info, ABR, sid, bufsize - ABR.span() * sizeof(value_type),
                                                     (void *)(w.data()));
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
      UnmanagedViewType<ordinal_type_array> ipiv(piv + rbeg * 4, 4 * m);

      const ordinal_type n = s.n - s.m;
      const ordinal_type nrhs = info.x.extent(1);
      const size_type bufsize_required = n * nrhs * sizeof(value_type);

      TACHO_TEST_FOR_ABORT(bufsize < bufsize_required, "bufsize is smaller than required");

      UnmanagedViewType<value_type_matrix> xB((value_type *)buf, n, nrhs);

      LDL_Supernodes<Algo::Workflow::Serial>::solve_lower(member, info, ipiv, xB, sid);

      CholSupernodes<Algo::Workflow::Serial>::update_solve_lower(member, info, xB, sid);
    }
    return 0;
  }

  template <typename MemberType, typename SupernodeInfoType>
  KOKKOS_INLINE_FUNCTION static int
  solve_upper_recursive_serial(MemberType &member, const SupernodeInfoType &info, const ordinal_type sid,
                               const bool final, typename SupernodeInfoType::ordinal_type_array::pointer_type piv,
                               typename SupernodeInfoType::value_type_array::pointer_type diag,
                               typename SupernodeInfoType::value_type_array::pointer_type buf,
                               const ordinal_type bufsize) {
    using supernode_info_type = SupernodeInfoType;
    using value_type = typename supernode_info_type::value_type;
    using value_type_matrix = typename supernode_info_type::value_type_matrix;
    using ordinal_type_array = typename supernode_info_type::ordinal_type_array;

    const auto &s = info.supernodes(sid);
    {
      const ordinal_type m = s.m;
      const ordinal_type rbeg = s.row_begin;
      UnmanagedViewType<ordinal_type_array> ipiv(piv + rbeg * 4, 4 * m);
      UnmanagedViewType<value_type_matrix> dblk(diag + rbeg * 2, m, 2);

      const ordinal_type n = s.n - s.m;
      const ordinal_type nrhs = info.x.extent(1);
      const ordinal_type bufsize_required = n * nrhs * sizeof(value_type);

      TACHO_TEST_FOR_ABORT(bufsize < bufsize_required, "bufsize is smaller than required");

      UnmanagedViewType<value_type_matrix> xB((value_type *)buf, n, nrhs);

      CholSupernodes<Algo::Workflow::Serial>::update_solve_upper(member, info, xB, sid);

      LDL_Supernodes<Algo::Workflow::Serial>::solve_upper(member, info, ipiv, dblk, xB, sid);
    }

    if (final) {
      // serial recursion
      for (ordinal_type i = 0; i < s.nchildren; ++i)
        solve_upper_recursive_serial(member, info, s.children[i], final, piv, diag, buf, bufsize);
    }
    return 0;
  }
};
} // namespace Tacho

#endif
