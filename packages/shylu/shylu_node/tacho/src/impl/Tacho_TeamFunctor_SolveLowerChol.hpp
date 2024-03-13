// clang-format off
/* =====================================================================================
Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.

SCR#:2790.0

This file is part of Tacho. Tacho is open source software: you can redistribute it
and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Kyungjoo Kim at <kyukim@sandia.gov,https://github.com/kyungjoo-kim>

Sandia National Laboratories, Albuquerque, NM, USA
===================================================================================== */
// clang-format on
#ifndef __TACHO_TEAMFUNCTOR_SOLVE_LOWER_CHOL_HPP__
#define __TACHO_TEAMFUNCTOR_SOLVE_LOWER_CHOL_HPP__

/// \file Tacho_TeamFunctor_SolveLowerChol.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

#include "Tacho_SupernodeInfo.hpp"

namespace Tacho {

template <typename SupernodeInfoType> struct TeamFunctor_SolveLowerChol {
public:
  typedef Kokkos::pair<ordinal_type, ordinal_type> range_type;

  typedef SupernodeInfoType supernode_info_type;
  typedef typename supernode_info_type::supernode_type supernode_type;
  typedef typename supernode_info_type::supernode_type_array supernode_type_array;

  typedef typename supernode_info_type::ordinal_type_array ordinal_type_array;
  typedef typename supernode_info_type::size_type_array size_type_array;

  typedef typename supernode_info_type::ordinal_pair_type_array ordinal_pair_type_array;

  typedef typename supernode_info_type::value_type value_type;
  typedef typename supernode_info_type::value_type_array value_type_array;
  typedef typename supernode_info_type::value_type_matrix value_type_matrix;

private:
  ConstUnmanagedViewType<supernode_type_array> _supernodes;
  ConstUnmanagedViewType<ordinal_pair_type_array> _sid_block_colidx;
  ConstUnmanagedViewType<ordinal_type_array> _gid_colidx;

  ConstUnmanagedViewType<ordinal_type_array> _compute_mode, _level_sids;
  ordinal_type _pbeg, _pend;

  UnmanagedViewType<value_type_matrix> _t;
  ordinal_type _nrhs;

  UnmanagedViewType<size_type_array> _buf_ptr;
  UnmanagedViewType<value_type_array> _buf;

public:
  KOKKOS_INLINE_FUNCTION
  TeamFunctor_SolveLowerChol() = delete;

  KOKKOS_INLINE_FUNCTION
  TeamFunctor_SolveLowerChol(const supernode_info_type &info, const ordinal_type_array &compute_mode,
                             const ordinal_type_array &level_sids, const value_type_matrix t,
                             const value_type_array buf)
      : _supernodes(info.supernodes), _sid_block_colidx(info.sid_block_colidx), _gid_colidx(info.gid_colidx),
        _compute_mode(compute_mode), _level_sids(level_sids), _t(t), _nrhs(t.extent(1)), _buf(buf) {}

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
    using GemvAlgoType = typename GemvAlgorithm_Team::type;
    using TrsvAlgoType = typename TrsvAlgorithm_Team::type;

    const value_type minus_one(-1), zero(0);
    {
      const ordinal_type m = s.m, n = s.n, n_m = n - m;
      if (m > 0) {
        value_type *aptr = s.u_buf;
        // solve
        UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
        aptr += m * m;

        const ordinal_type offm = s.row_begin;
        auto tT = Kokkos::subview(_t, range_type(offm, offm + m), Kokkos::ALL());
        Trsv<Uplo::Upper, Trans::ConjTranspose, TrsvAlgoType>::invoke(member, Diag::NonUnit(), ATL, tT);

        if (n_m > 0) {
          // update
          member.team_barrier();
          UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n;
          UnmanagedViewType<value_type_matrix> bB(bptr, n_m, _nrhs);
          Gemv<Trans::ConjTranspose, GemvAlgoType>::invoke(member, minus_one, ATR, tT, zero, bB);
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

    const value_type minus_one(-1), one(1), zero(0);
    {
      const ordinal_type m = s.m, n = s.n, n_m = n - m;
      if (m > 0) {
        value_type *aptr = s.u_buf;

        UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
        aptr += m * m;
        UnmanagedViewType<value_type_matrix> bT(bptr, m, _nrhs);
        bptr += m * _nrhs;

        const ordinal_type offm = s.row_begin;
        auto tT = Kokkos::subview(_t, range_type(offm, offm + m), Kokkos::ALL());

        Gemv<Trans::ConjTranspose, GemvAlgoType>::invoke(member, one, ATL, tT, zero, bT);

        if (n_m > 0) {
          // solve offdiag
          member.team_barrier();
          UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m);
          UnmanagedViewType<value_type_matrix> bB(bptr, n_m, _nrhs);

          Gemv<Trans::ConjTranspose, GemvAlgoType>::invoke(member, minus_one, ATR, bT, zero, bB);
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

        // copy to t
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(member, m * _nrhs),
            [&, m](const ordinal_type &k) { // Value capture is a workaround for cuda + gcc-7.2 compiler bug w/c++14
              const ordinal_type i = k % m, j = k / m;
              tT(i, j) = bT(i, j);
            });

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
      if (m > 0 && n > 0) {
        value_type *aptr = s.u_buf;

        UnmanagedViewType<value_type_matrix> AT(aptr, m, n);
        UnmanagedViewType<value_type_matrix> b(bptr, n, _nrhs);

        const ordinal_type offm = s.row_begin;
        auto tT = Kokkos::subview(_t, range_type(offm, offm + m), Kokkos::ALL());

        Gemv<Trans::ConjTranspose, GemvAlgoType>::invoke(member, one, AT, tT, zero, b);
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
      typedef SolveTag<Var> solve_tag_type;
      const supernode_type &s = _supernodes(sid);
      value_type *bptr = _buf.data() + _buf_ptr(member.league_rank());
      if (solve_tag_type::variant == 0)
        solve_var0(member, s, bptr);
      else if (solve_tag_type::variant == 1)
        solve_var1(member, s, bptr);
      else if (solve_tag_type::variant == 2)
        solve_var2(member, s, bptr);
      else
        Kokkos::printf("Error: TeamFunctorSolveLowerChol::SolveTag, algorithm variant is not supported\n");
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
      const ordinal_type sid = _level_sids(p);
      typedef UpdateTag<Var> update_tag_type;
      const supernode_type &s = _supernodes(sid);
      value_type *bptr = _buf.data() + _buf_ptr(member.league_rank());
      if (update_tag_type::variant == 0)
        update_var0(member, s, bptr);
      else if (update_tag_type::variant == 1)
        update_var1(member, s, bptr);
      else if (update_tag_type::variant == 2)
        update_var2(member, s, bptr);
      else
        Kokkos::printf("Error: TeamFunctorSolveLowerChol::UpdateTag, algorithm variant is not supported\n");
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
