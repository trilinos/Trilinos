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
#ifndef __TACHO_TEAMFUNCTOR_SOLVE_UPPER_LDL_HPP__
#define __TACHO_TEAMFUNCTOR_SOLVE_UPPER_LDL_HPP__

/// \file Tacho_TeamFunctor_SolveUpperLDL.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

#include "Tacho_SupernodeInfo.hpp"

namespace Tacho {

template <typename SupernodeInfoType> struct TeamFunctor_SolveUpperLDL {
public:
  using range_type = Kokkos::pair<ordinal_type, ordinal_type>;

  using supernode_info_type = SupernodeInfoType;
  using supernode_type = typename supernode_info_type::supernode_type;
  using supernode_type_array = typename supernode_info_type::supernode_type_array;

  using ordinal_type_array = typename supernode_info_type::ordinal_type_array;
  using size_type_array = typename supernode_info_type::size_type_array;

  using ordinal_pair_type_array = typename supernode_info_type::ordinal_pair_type_array;

  using value_type = typename supernode_info_type::value_type;
  using value_type_array = typename supernode_info_type::value_type_array;
  using value_type_matrix = typename supernode_info_type::value_type_matrix;

private:
  ConstUnmanagedViewType<supernode_type_array> _supernodes;
  ConstUnmanagedViewType<ordinal_type_array> _gid_colidx;

  ConstUnmanagedViewType<ordinal_type_array> _compute_mode, _level_sids;
  ordinal_type _pbeg, _pend;

  ConstUnmanagedViewType<ordinal_type_array> _piv;
  ConstUnmanagedViewType<value_type_array> _diag;
  UnmanagedViewType<value_type_matrix> _t;
  ordinal_type _nrhs;

  UnmanagedViewType<size_type_array> _buf_ptr;
  UnmanagedViewType<value_type_array> _buf;

public:
  KOKKOS_INLINE_FUNCTION
  TeamFunctor_SolveUpperLDL() = delete;

  KOKKOS_INLINE_FUNCTION
  TeamFunctor_SolveUpperLDL(const supernode_info_type &info, const ordinal_type_array &compute_mode,
                            const ordinal_type_array &level_sids, const ordinal_type_array &piv,
                            const value_type_array &diag, const value_type_matrix &t, const value_type_array &buf)
      : _supernodes(info.supernodes), _gid_colidx(info.gid_colidx), _compute_mode(compute_mode),
        _level_sids(level_sids), _piv(piv), _diag(diag), _t(t), _nrhs(t.extent(1)), _buf(buf) {}

  inline void setRange(const ordinal_type pbeg, const ordinal_type pend) {
    _pbeg = pbeg;
    _pend = pend;
  }

  inline void setBufferPtr(const size_type_array &buf_ptr) { _buf_ptr = buf_ptr; }

  ///
  /// Algorithm Variant 0: gemv - trsv
  ///
  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void solve_var0(MemberType &member, const supernode_type &s, value_type *bptr) const {
    using TrsvAlgoType = typename TrsvAlgorithm::type;
    using GemvAlgoType = typename GemvAlgorithm::type;

    const value_type minus_one(-1), one(1);
    {
      const ordinal_type m = s.m, n = s.n, n_m = n - m;
      if (m > 0) {
        value_type *aptr = s.u_buf;
        // solve
        const UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
        aptr += m * m;
        const ordinal_type offm = s.row_begin;
        const auto tT = Kokkos::subview(_t, range_type(offm, offm + m), Kokkos::ALL());

        ConstUnmanagedViewType<ordinal_type_array> P(_piv.data() + offm * 4, m * 4);
        ConstUnmanagedViewType<value_type_matrix> D(_diag.data() + offm * 2, m, 2);

        Scale2x2_BlockInverseDiagonals<Side::Left, Algo::Internal>::invoke(member, P, D, tT);

        if (n_m > 0) {
          // update
          const UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n;
          const UnmanagedViewType<value_type_matrix> bB(bptr, n_m, _nrhs);
          Gemv<Trans::NoTranspose, GemvAlgoType>::invoke(member, minus_one, ATR, bB, one, tT);
          member.team_barrier();
        }
        Trsv<Uplo::Lower, Trans::Transpose, TrsvAlgoType>::invoke(member, Diag::Unit(), ATL, tT);

        if (!s.do_not_apply_pivots) {
          ConstUnmanagedViewType<ordinal_type_array> fpiv(P.data() + m, m);
          ApplyPivots<PivotMode::Flame, Side::Left, Direct::Backward, Algo::Internal> /// row inter-change
              ::invoke(member, fpiv, tT);
        }
      }
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void update_var0(MemberType &member, const supernode_type &s, value_type *bptr) const {
    {
      const ordinal_type m = s.m, n = s.n, n_m = n - m;
      if (n_m > 0) {
        // update
        const UnmanagedViewType<value_type_matrix> bB(bptr, n_m, _nrhs);
        const ordinal_type goffset = s.gid_col_begin + s.m;
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(member, n_m),
            [&,
             goffset](const ordinal_type &i) { // Value capture is a workaround for cuda + gcc-7.2 compiler bug w/c++14
              // for (ordinal_type i=0;i<n;++i) {
              const ordinal_type row = _gid_colidx(i + goffset);
              for (ordinal_type j = 0; j < _nrhs; ++j)
                bB(i, j) = _t(row, j);
            });
      }
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void solve_var1(MemberType &member, const supernode_type &s, value_type *bptr) const {
    using GemvAlgoType = typename GemvAlgorithm::type;

    const value_type minus_one(-1), one(1), zero(0);
    {
      const ordinal_type m = s.m, n = s.n, n_m = n - m;
      if (m > 0) {
        value_type *aptr = s.u_buf;
        // solve
        const UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
        aptr += ATL.span();
        const UnmanagedViewType<value_type_matrix> bT(bptr, m, _nrhs);
        bptr += bT.span();

        const ordinal_type offm = s.row_begin;
        const auto tT = Kokkos::subview(_t, range_type(offm, offm + m), Kokkos::ALL());

        ConstUnmanagedViewType<ordinal_type_array> P(_piv.data() + offm * 4, m * 4);
        ConstUnmanagedViewType<value_type_matrix> D(_diag.data() + offm * 2, m, 2);

        Scale2x2_BlockInverseDiagonals<Side::Left, Algo::Internal>::invoke(member, P, D, tT);

        if (n_m > 0) {
          // update
          const UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n;
          const UnmanagedViewType<value_type_matrix> bB(bptr, n_m, _nrhs);
          Gemv<Trans::NoTranspose, GemvAlgoType>::invoke(member, minus_one, ATR, bB, one, tT);
          member.team_barrier();
        }
        Gemv<Trans::Transpose, GemvAlgoType>::invoke(member, one, ATL, tT, zero, bT);
        member.team_barrier();

        if (s.do_not_apply_pivots) {
          Copy<Algo::Internal>::invoke(member, tT, bT);
        } else {
          ConstUnmanagedViewType<ordinal_type_array> peri(P.data() + 3 * m, m);
          ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::Internal>::invoke(member, bT, peri, tT);
        }
      }
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void update_var1(MemberType &member, const supernode_type &s, value_type *bptr) const {
    {
      const ordinal_type m = s.m, n = s.n, n_m = n - m;
      if (n_m > 0) {
        UnmanagedViewType<value_type_matrix> bB(bptr + m * _nrhs, n_m, _nrhs);

        // update
        const ordinal_type goffset = s.gid_col_begin + s.m;
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(member, n_m),
            [&,
             goffset](const ordinal_type &i) { // Value capture is a workaround for cuda + gcc-7.2 compiler bug w/c++14
              // for (ordinal_type i=0;i<n;++i) {
              const ordinal_type row = _gid_colidx(i + goffset);
              for (ordinal_type j = 0; j < _nrhs; ++j)
                bB(i, j) = _t(row, j);
            });
      }
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void solve_var2(MemberType &member, const supernode_type &s, value_type *bptr) const {
    using GemvAlgoType = typename GemvAlgorithm::type;

    const ordinal_type nrhs = _t.extent(1);
    const value_type one(1), zero(0);
    {
      const ordinal_type m = s.m, n = s.n;
      if (m > 0 && n > 0) {
        value_type *aptr = s.u_buf;

        const UnmanagedViewType<value_type_matrix> AT(aptr, m, n);
        const UnmanagedViewType<value_type_matrix> b(bptr, n, nrhs);

        const ordinal_type offm = s.row_begin;
        const auto tT = Kokkos::subview(_t, range_type(offm, offm + m), Kokkos::ALL());
        const UnmanagedViewType<value_type_matrix> bT(bptr, m, nrhs);

        ConstUnmanagedViewType<ordinal_type_array> P(_piv.data() + offm * 4, m * 4);
        ConstUnmanagedViewType<value_type_matrix> D(_diag.data() + offm * 2, m, 2);

        Scale2x2_BlockInverseDiagonals<Side::Left, Algo::Internal>::invoke(member, P, D, bT);
        member.team_barrier();

        Gemv<Trans::NoTranspose, GemvAlgoType>::invoke(member, one, AT, b, zero, tT);
        member.team_barrier();

        if (!s.do_not_apply_pivots) {
          Copy<Algo::Internal>::invoke(member, bT, tT);
          member.team_barrier();

          ConstUnmanagedViewType<ordinal_type_array> peri(P.data() + 3 * m, m);
          ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::Internal>::invoke(member, bT, peri, tT);
        }
      }
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void update_var2(MemberType &member, const supernode_type &s, value_type *bptr) const {
    {
      const ordinal_type m = s.m, n = s.n;
      UnmanagedViewType<value_type_matrix> b(bptr, n, _nrhs);
      if (n > 0) {
        const ordinal_type offm = s.row_begin;
        const ordinal_type goffset = s.gid_col_begin + s.m;
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(member, n),
            [&, m, goffset,
             offm](const ordinal_type &i) { // Value capture is a workaround for cuda + gcc-7.2 compiler bug w/c++14
              for (ordinal_type j = 0; j < _nrhs; ++j) {
                if (i < m) {
                  b(i, j) = _t(offm + i, j);
                } else {
                  const ordinal_type row = _gid_colidx(i - m + goffset);
                  b(i, j) = _t(row, j);
                }
              }
            });
      }
    }
  }

  template <int Var> struct SolveTag {
    enum { variant = Var };
  };
  template <int Var> struct UpdateTag {
    enum { variant = Var };
  };
  struct DummyTag {};

  template <typename MemberType, int Var>
  KOKKOS_INLINE_FUNCTION void operator()(const SolveTag<Var> &, const MemberType &member) const {
    const ordinal_type p = _pbeg + member.league_rank();
    const ordinal_type sid = _level_sids(p);
    const ordinal_type mode = _compute_mode(sid);
    if (p < _pend && mode == 1) {
      using solve_tag_type = SolveTag<Var>;

      const supernode_type &s = _supernodes(sid);
      value_type *bptr = _buf.data() + _buf_ptr(member.league_rank());
      if (solve_tag_type::variant == 0) {
        solve_var0(member, s, bptr);
      } else if (solve_tag_type::variant == 1) {
        solve_var1(member, s, bptr);
      } else if (solve_tag_type::variant == 2 || solve_tag_type::variant == 3) {
        solve_var2(member, s, bptr);
      }
    } else if (mode == -1) {
      Kokkos::printf("Error: TeamFunctorSolveUpperChol::SolveTag, computing mode is not determined\n");
    } else {
      // skip
    }
  }

  template <typename MemberType, int Var>
  KOKKOS_INLINE_FUNCTION void operator()(const UpdateTag<Var> &, const MemberType &member) const {
    const ordinal_type p = _pbeg + member.league_rank();
    if (p < _pend) {
      using update_tag_type = UpdateTag<Var>;

      const ordinal_type sid = _level_sids(p);
      const supernode_type &s = _supernodes(sid);
      value_type *bptr = _buf.data() + _buf_ptr(member.league_rank());

      if (update_tag_type::variant == 0) {
        update_var0(member, s, bptr);
      } else if (update_tag_type::variant == 1) {
        update_var1(member, s, bptr);
      } else if (update_tag_type::variant == 2 || update_tag_type::variant == 3) {
        update_var2(member, s, bptr);
      }
    } else {
      // skip
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const DummyTag &, const MemberType &member) const {}
};

} // namespace Tacho

#endif
