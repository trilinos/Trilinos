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
#ifndef __TACHO_TEAMFUNCTOR_SOLVE_LOWER_LU_HPP__
#define __TACHO_TEAMFUNCTOR_SOLVE_LOWER_LU_HPP__

/// \file Tacho_TeamFunctor_SolveLowerLU.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

#include "Tacho_SupernodeInfo.hpp"

namespace Tacho {

template <typename SupernodeInfoType> struct TeamFunctor_SolveLowerLU {
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
  ConstUnmanagedViewType<ordinal_pair_type_array> _sid_block_colidx;
  ConstUnmanagedViewType<ordinal_type_array> _gid_colidx;

  ConstUnmanagedViewType<ordinal_type_array> _compute_mode, _level_sids;
  ordinal_type _pbeg, _pend;

  ConstUnmanagedViewType<ordinal_type_array> _piv;

  UnmanagedViewType<value_type_matrix> _t;
  ordinal_type _nrhs;

  UnmanagedViewType<size_type_array> _buf_ptr;
  UnmanagedViewType<value_type_array> _buf;

public:
  KOKKOS_INLINE_FUNCTION
  TeamFunctor_SolveLowerLU() = delete;

  KOKKOS_INLINE_FUNCTION
  TeamFunctor_SolveLowerLU(const supernode_info_type &info, const ordinal_type_array &compute_mode,
                           const ordinal_type_array &level_sids, const ordinal_type_array &piv,
                           const value_type_matrix t, const value_type_array buf)
      : _supernodes(info.supernodes), _sid_block_colidx(info.sid_block_colidx), _gid_colidx(info.gid_colidx),
        _compute_mode(compute_mode), _level_sids(level_sids), _piv(piv), _t(t), _nrhs(t.extent(1)), _buf(buf) {}

  inline void setRange(const ordinal_type pbeg, const ordinal_type pend) {
    _pbeg = pbeg;
    _pend = pend;
  }

  inline void setBufferPtr(const size_type_array &buf_ptr) { _buf_ptr = buf_ptr; }

  ///
  /// Algorithm Variant 0: trsv - gemv
  ///
  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void solve_var0(MemberType &member, const supernode_type &s, value_type *bptr) const {
    using TrsvAlgoType = typename TrsvAlgorithm_Team::type;
    using GemvAlgoType = typename GemvAlgorithm_Team::type;

    const value_type minus_one(-1), zero(0);
    {
      const ordinal_type m = s.m, n = s.n, n_m = n - m;
      if (m > 0) {
        // solve
        UnmanagedViewType<value_type_matrix> ATL(s.u_buf, m, m);

        const ordinal_type offm = s.row_begin;
        auto tT = Kokkos::subview(_t, range_type(offm, offm + m), Kokkos::ALL());
        auto fpiv = ConstUnmanagedViewType<ordinal_type_array>(_piv.data() + 4 * offm + m, m);

        if (!s.do_not_apply_pivots) {
          ApplyPivots<PivotMode::Flame, Side::Left, Direct::Forward, Algo::Internal> /// row inter-change
              ::invoke(member, fpiv, tT);
        }
        Trsv<Uplo::Lower, Trans::NoTranspose, TrsvAlgoType>::invoke(member, Diag::Unit(), ATL, tT);

        if (n_m > 0) {
          // update
          member.team_barrier();
          UnmanagedViewType<value_type_matrix> AL(s.l_buf, n, m);
          const auto ABL = Kokkos::subview(AL, range_type(m, n), Kokkos::ALL());
          UnmanagedViewType<value_type_matrix> bB(bptr, n_m, _nrhs);
          Gemv<Trans::NoTranspose, GemvAlgoType>::invoke(member, minus_one, ABL, tT, zero, bB);
        }
      }
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void update_var0(MemberType &member, const supernode_type &s, value_type *bptr) const {
    {
      const ordinal_type m = s.m, n = s.n, n_m = n - m;
      if (n_m > 0) {
        UnmanagedViewType<value_type_matrix> bB(bptr, n_m, _nrhs);
        // update
        const ordinal_type sbeg = s.sid_col_begin + 1, send = s.sid_col_end - 1;
        for (ordinal_type i = sbeg, ip = 0 /*is=0*/; i < send; ++i) {
          const ordinal_type tbeg = _sid_block_colidx(i).second, tend = _sid_block_colidx(i + 1).second,
                             tcnt = tend - tbeg;

          Kokkos::parallel_for(
              Kokkos::TeamVectorRange(member, tcnt),
              [&, tbeg,
               ip](const ordinal_type &ii) { // Value capture is a workaround for cuda + gcc-7.2 compiler bug w/c++14
                const ordinal_type it = tbeg + ii;
                const ordinal_type is = ip + ii;
                // for (ordinal_type it=tbeg;it<tend;++it,++is) {
                const ordinal_type row = _gid_colidx(s.gid_col_begin + it);
                for (ordinal_type j = 0; j < _nrhs; ++j)
                  Kokkos::atomic_add(&_t(row, j), bB(is, j));
              });
          ip += tcnt;
        }
      }
    }
  }

  ///
  /// Algorithm Variant 1: gemv - gemv
  ///
  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void solve_var1(MemberType &member, const supernode_type &s, value_type *bptr) const {
    using GemvAlgoType = typename GemvAlgorithm::type;

    const value_type one(1), minus_one(-1), zero(0);
    {
      const ordinal_type m = s.m, n = s.n, n_m = n - m;
      if (m > 0) {
        // solve
        UnmanagedViewType<value_type_matrix> AL(s.l_buf, n, m);
        const auto ATL = Kokkos::subview(AL, range_type(0, m), Kokkos::ALL());
        UnmanagedViewType<value_type_matrix> bT(bptr, m, _nrhs);

        const ordinal_type offm = s.row_begin;
        const auto tT = Kokkos::subview(_t, range_type(offm, offm + m), Kokkos::ALL());
        const ConstUnmanagedViewType<ordinal_type_array> perm(_piv.data() + 4 * offm + 2 * m, m);

        if (s.do_not_apply_pivots) {
          Copy<Algo::Internal>::invoke(member, bT, tT);
        } else {
          ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::Internal>::invoke(member, tT, perm, bT);
        }
        member.team_barrier();

        Gemv<Trans::NoTranspose, GemvAlgoType>::invoke(member, one, ATL, bT, zero, tT);

        if (n_m > 0) {
          // update
          member.team_barrier();
          const auto ABL = Kokkos::subview(AL, range_type(m, n), Kokkos::ALL());
          UnmanagedViewType<value_type_matrix> bB(bptr + bT.span(), n_m, _nrhs);
          Gemv<Trans::NoTranspose, GemvAlgoType>::invoke(member, minus_one, ABL, tT, zero, bB);
        }
      }
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void update_var1(MemberType &member, const supernode_type &s, value_type *bptr) const {
    {
      const ordinal_type m = s.m, n = s.n, n_m = n - m;
      if (m > 0) {
        UnmanagedViewType<value_type_matrix> bT(bptr, m, _nrhs);
        bptr += m * _nrhs;

        const ordinal_type offm = s.row_begin;
        auto tT = Kokkos::subview(_t, range_type(offm, offm + m), Kokkos::ALL());

        if (n_m > 0) {
          UnmanagedViewType<value_type_matrix> bB(bptr, n_m, _nrhs);

          // update
          const ordinal_type sbeg = s.sid_col_begin + 1, send = s.sid_col_end - 1;
          for (ordinal_type i = sbeg, ip = 0 /*is=0*/; i < send; ++i) {
            const ordinal_type tbeg = _sid_block_colidx(i).second, tend = _sid_block_colidx(i + 1).second,
                               tcnt = tend - tbeg;

            Kokkos::parallel_for(
                Kokkos::TeamVectorRange(member, tcnt),
                [&, tbeg](
                    const ordinal_type &ii) { // Value capture is a workaround for cuda + gcc-7.2 compiler bug w/c++14
                  const ordinal_type it = tbeg + ii;
                  const ordinal_type is = ip + ii;
                  // for (ordinal_type it=tbeg;it<tend;++it,++is) {
                  const ordinal_type row = _gid_colidx(s.gid_col_begin + it);
                  for (ordinal_type j = 0; j < _nrhs; ++j)
                    Kokkos::atomic_add(&_t(row, j), bB(is, j));
                });
            ip += tcnt;
          }
        }
      }
    }
  }

  ///
  /// Algorithm Variant 2: gemv
  ///
  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void solve_var2(MemberType &member, const supernode_type &s, value_type *bptr) const {
    using GemvAlgoType = typename GemvAlgorithm::type;

    const value_type one(1), zero(0);
    {
      const ordinal_type m = s.m, n = s.n;
      if (m > 0) {
        // solve
        UnmanagedViewType<value_type_matrix> AL(s.l_buf, n, m);
        UnmanagedViewType<value_type_matrix> b(bptr, n, _nrhs);

        const ordinal_type offm = s.row_begin;
        const auto tT = Kokkos::subview(_t, range_type(offm, offm + m), Kokkos::ALL());

        if (!s.do_not_apply_pivots) {
          UnmanagedViewType<value_type_matrix> bT(bptr, m, _nrhs);
          ConstUnmanagedViewType<ordinal_type_array> perm(_piv.data() + 4 * offm + 2 * m, m);
          Copy<Algo::Internal>::invoke(member, bT, tT);
          member.team_barrier();
          ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::Internal>::invoke(member, bT, perm, tT);
          member.team_barrier();
        }
        Gemv<Trans::NoTranspose, GemvAlgoType>::invoke(member, one, AL, tT, zero, b);
      }
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void update_var2(MemberType &member, const supernode_type &s, value_type *bptr) const {
    {
      const ordinal_type m = s.m, n = s.n, n_m = n - m;
      UnmanagedViewType<value_type_matrix> b(bptr, n, _nrhs);
      if (m > 0) {
        const ordinal_type offm = s.row_begin;
        auto tT = Kokkos::subview(_t, range_type(offm, offm + m), Kokkos::ALL());

        // copy to t
        Kokkos::parallel_for(Kokkos::TeamVectorRange(member, m * _nrhs),
                             [&, m](const ordinal_type &k) { /// compiler bug with c++14 lambda capturing and workaround
                               const ordinal_type i = k % m, j = k / m;
                               tT(i, j) = b(i, j);
                             });
        member.team_barrier();

        if (n_m > 0) {
          // update
          const ordinal_type sbeg = s.sid_col_begin + 1, send = s.sid_col_end - 1;
          for (ordinal_type i = sbeg, ip = 0 /*is=0*/; i < send; ++i) {
            const ordinal_type tbeg = _sid_block_colidx(i).second, tend = _sid_block_colidx(i + 1).second,
                               tcnt = tend - tbeg;

            Kokkos::parallel_for(
                Kokkos::TeamVectorRange(member, tcnt),
                [&, ip, m, tbeg](const ordinal_type &ii) { /// compiler bug with c++14 lambda capturing and workaround
                  const ordinal_type it = tbeg + ii;
                  const ordinal_type is = ip + ii;
                  // for (ordinal_type it=tbeg;it<tend;++it,++is) {
                  const ordinal_type row = _gid_colidx(s.gid_col_begin + it);
                  for (ordinal_type j = 0; j < _nrhs; ++j)
                    Kokkos::atomic_add(&_t(row, j), b(is + m, j));
                });
            ip += tcnt;
          }
        }
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
    }
    if (mode == -1) {
      Kokkos::printf("Error: TeamFunctorSolveLowerChol::SolveTag, computing mode is not determined\n");
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
  KOKKOS_INLINE_FUNCTION void operator()(const DummyTag &, const MemberType &member) const {
    // do nothing
  }
};
} // namespace Tacho

#endif
