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
#ifndef __TACHO_TEAMFUNCTOR_FACTORIZE_LU_HPP__
#define __TACHO_TEAMFUNCTOR_FACTORIZE_LU_HPP__

/// \file Tacho_TeamFunctor_FactorizeLU.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

#include "Tacho_SupernodeInfo.hpp"

namespace Tacho {

template <typename SupernodeInfoType> struct TeamFunctor_FactorizeLU {
public:
  using range_type = Kokkos::pair<ordinal_type, ordinal_type>;

  using supernode_info_type = SupernodeInfoType;
  using supernode_type = typename supernode_info_type::supernode_type;

  using ordinal_type_array = typename supernode_info_type::ordinal_type_array;
  using size_type_array = typename supernode_info_type::size_type_array;

  using value_type = typename supernode_info_type::value_type;
  using value_type_array = typename supernode_info_type::value_type_array;
  using value_type_matrix = typename supernode_info_type::value_type_matrix;

private:
  supernode_info_type _info;
  ordinal_type_array _compute_mode, _level_sids;
  ordinal_type _pbeg, _pend;

  ordinal_type_array _piv;

  size_type_array _buf_ptr;
  value_type_array _buf;

  int *_rval;

public:
  KOKKOS_INLINE_FUNCTION
  TeamFunctor_FactorizeLU() = delete;

  KOKKOS_INLINE_FUNCTION
  TeamFunctor_FactorizeLU(const supernode_info_type &info, const ordinal_type_array &compute_mode,
                          const ordinal_type_array &level_sids, const ordinal_type_array &piv,
                          const value_type_array buf, int *rval)
      : _info(info), _compute_mode(compute_mode), _level_sids(level_sids), _piv(piv), _buf(buf), _rval(rval) {}

  inline void setRange(const ordinal_type pbeg, const ordinal_type pend) {
    _pbeg = pbeg;
    _pend = pend;
  }

  inline void setBufferPtr(const size_type_array &buf_ptr) { _buf_ptr = buf_ptr; }

  ///
  /// Main functions
  ///
  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void factorize_var0(MemberType &member, const supernode_type &s, const ordinal_type_array &P,
                                             const value_type_matrix &ABR) const {
    using LU_AlgoType = typename LU_Algorithm_Team::type;
    using GemmAlgoType = typename GemmAlgorithm_Team::type;
    using TrsmAlgoType = typename TrsmAlgorithm_Team::type;

    int err = 0;
    const ordinal_type m = s.m, n = s.n, n_m = n - m;
    if (m > 0) {
      UnmanagedViewType<value_type_matrix> AT(s.u_buf, m, n);

      err = LU<LU_AlgoType>::invoke(member, AT, P);
      member.team_barrier();
      if (err != 0) {
        Kokkos::atomic_add(_rval, 1);
        return;
      }

      LU<LU_AlgoType>::modify(member, m, P);
      member.team_barrier();

      if (n_m > 0) {
        const value_type one(1), minus_one(-1), zero(0);
        UnmanagedViewType<value_type_matrix> ATL(s.u_buf, m, m);
        UnmanagedViewType<value_type_matrix> ATR(s.u_buf + ATL.span(), m, n_m);
        UnmanagedViewType<value_type_matrix> AL(s.l_buf, n, m);
        const auto ABL = Kokkos::subview(AL, range_type(m, n), Kokkos::ALL());

        Trsm<Side::Right, Uplo::Upper, Trans::NoTranspose, TrsmAlgoType>::invoke(member, Diag::NonUnit(), one, ATL,
                                                                                 ABL);
        member.team_barrier();

        Gemm<Trans::NoTranspose, Trans::NoTranspose, GemmAlgoType>::invoke(member, minus_one, ABL, ATR, zero, ABR);
      }
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void factorize_var1(MemberType &member, const supernode_type &s, const ordinal_type_array &P,
                                             const value_type_matrix &T, const value_type_matrix &ABR) const {
    using LU_AlgoType = typename LU_Algorithm::type;
    using TrsmAlgoType = typename TrsmAlgorithm::type;
    using GemmAlgoType = typename GemmAlgorithm::type;

    int err = 0;
    const value_type one(1), minus_one(-1), zero(0);
    const ordinal_type m = s.m, n = s.n, n_m = n - m;
    if (m > 0) {
      UnmanagedViewType<value_type_matrix> AT(s.u_buf, m, n);

      err = LU<LU_AlgoType>::invoke(member, AT, P);
      member.team_barrier();
      if (err != 0) {
        Kokkos::atomic_add(_rval, 1);
        return;
      }

      LU<LU_AlgoType>::modify(member, m, P);
      member.team_barrier();

      if (n_m > 0) {
        UnmanagedViewType<value_type_matrix> ATL(s.u_buf, m, m);
        UnmanagedViewType<value_type_matrix> ATR(s.u_buf + ATL.span(), m, n_m);
        UnmanagedViewType<value_type_matrix> AL(s.l_buf, n, m);
        const auto ATL2 = Kokkos::subview(AL, range_type(0, m), Kokkos::ALL());
        const auto ABL = Kokkos::subview(AL, range_type(m, n), Kokkos::ALL());

        Copy<Algo::Internal>::invoke(member, T, ATL);
        Trsm<Side::Right, Uplo::Upper, Trans::NoTranspose, TrsmAlgoType>::invoke(member, Diag::NonUnit(), one, ATL,
                                                                                 ABL);
        member.team_barrier();

        SetIdentity<Algo::Internal>::invoke(member, ATL, one);
        SetIdentity<Algo::Internal>::invoke(member, ATL2, one);
        member.team_barrier();

        Trsm<Side::Left, Uplo::Upper, Trans::NoTranspose, TrsmAlgoType>::invoke(member, Diag::NonUnit(), one, T, ATL);
        Trsm<Side::Left, Uplo::Lower, Trans::NoTranspose, TrsmAlgoType>::invoke(member, Diag::Unit(), one, T, ATL2);
        member.team_barrier();

        Gemm<Trans::NoTranspose, Trans::NoTranspose, GemmAlgoType>::invoke(member, minus_one, ABL, ATR, zero, ABR);
      } else {
        UnmanagedViewType<value_type_matrix> ATL(s.u_buf, m, m);
        UnmanagedViewType<value_type_matrix> ATL2(s.l_buf, m, m);

        Copy<Algo::Internal>::invoke(member, T, ATL);
        member.team_barrier();

        SetIdentity<Algo::Internal>::invoke(member, ATL, one);
        SetIdentity<Algo::Internal>::invoke(member, ATL2, one);
        member.team_barrier();

        Trsm<Side::Left, Uplo::Upper, Trans::NoTranspose, TrsmAlgoType>::invoke(member, Diag::NonUnit(), one, T, ATL);
        Trsm<Side::Left, Uplo::Lower, Trans::NoTranspose, TrsmAlgoType>::invoke(member, Diag::Unit(), one, T, ATL2);
      }
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void factorize_var2(MemberType &member, const supernode_type &s, const ordinal_type_array &P,
                                             const value_type_matrix &T, const value_type_matrix &ABR) const {
    using LU_AlgoType = typename LU_Algorithm::type;
    using TrsmAlgoType = typename TrsmAlgorithm::type;
    using GemmAlgoType = typename GemmAlgorithm::type;

    int err = 0;
    const value_type one(1), minus_one(-1), zero(0);
    const ordinal_type m = s.m, n = s.n, n_m = n - m;
    if (m > 0) {
      UnmanagedViewType<value_type_matrix> AT(s.u_buf, m, n);

      err = LU<LU_AlgoType>::invoke(member, AT, P);
      member.team_barrier();
      if (err != 0) {
        Kokkos::atomic_add(_rval, 1);
        return;
      }

      LU<LU_AlgoType>::modify(member, m, P);
      member.team_barrier();

      if (n_m > 0) {
        UnmanagedViewType<value_type_matrix> ATL(s.u_buf, m, m);
        UnmanagedViewType<value_type_matrix> ATR(s.u_buf + ATL.span(), m, n_m);
        UnmanagedViewType<value_type_matrix> AL(s.l_buf, n, m);
        const auto ATL2 = Kokkos::subview(AL, range_type(0, m), Kokkos::ALL());
        const auto ABL = Kokkos::subview(AL, range_type(m, n), Kokkos::ALL());

        Copy<Algo::Internal>::invoke(member, T, ATL);
        Trsm<Side::Right, Uplo::Upper, Trans::NoTranspose, TrsmAlgoType>::invoke(member, Diag::NonUnit(), one, ATL,
                                                                                 ABL);
        member.team_barrier();

        SetIdentity<Algo::Internal>::invoke(member, ATL, minus_one);
        SetIdentity<Algo::Internal>::invoke(member, ATL2, minus_one);
        member.team_barrier();

        Gemm<Trans::NoTranspose, Trans::NoTranspose, GemmAlgoType>::invoke(member, minus_one, ABL, ATR, zero, ABR);
        member.team_barrier();

        Trsm<Side::Left, Uplo::Upper, Trans::NoTranspose, TrsmAlgoType>::invoke(member, Diag::NonUnit(), minus_one, T,
                                                                                AT);
        Trsm<Side::Right, Uplo::Lower, Trans::NoTranspose, TrsmAlgoType>::invoke(member, Diag::Unit(), minus_one, T,
                                                                                 AL);
      } else {
        UnmanagedViewType<value_type_matrix> ATL(s.u_buf, m, m);
        UnmanagedViewType<value_type_matrix> ATL2(s.l_buf, m, m);

        Copy<Algo::Internal>::invoke(member, T, ATL);
        member.team_barrier();

        SetIdentity<Algo::Internal>::invoke(member, ATL, one);
        SetIdentity<Algo::Internal>::invoke(member, ATL2, one);
        member.team_barrier();

        Trsm<Side::Left, Uplo::Upper, Trans::NoTranspose, TrsmAlgoType>::invoke(member, Diag::NonUnit(), one, T, ATL);
        Trsm<Side::Left, Uplo::Lower, Trans::NoTranspose, TrsmAlgoType>::invoke(member, Diag::Unit(), one, T, ATL2);
      }
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void update(MemberType &member, const supernode_type &cur,
                                     const value_type_matrix &ABR) const {

    static constexpr bool runOnHost = run_tacho_on_host_v<typename value_type_matrix::execution_space>;

    const auto info = _info;
    value_type *buf = ABR.data() + ABR.span();
    const ordinal_type sbeg = cur.sid_col_begin + 1, send = cur.sid_col_end - 1;

    const ordinal_type srcbeg = info.sid_block_colidx(sbeg).second, srcend = info.sid_block_colidx(send).second,
                       srcsize = srcend - srcbeg;

    // short cut to direct update
    if ((send - sbeg) == 1) {
      const auto &s = info.supernodes(info.sid_block_colidx(sbeg).first);
      const ordinal_type tgtbeg = info.sid_block_colidx(s.sid_col_begin).second,
                         tgtend = info.sid_block_colidx(s.sid_col_end - 1).second, tgtsize = tgtend - tgtbeg;

      if (srcsize == tgtsize) {
        /* */ value_type *tgt = s.u_buf;
        const value_type *src = (value_type *)ABR.data();

        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(member, srcsize),
            [&, srcsize, src,
             tgt](const ordinal_type &j) { // Value capture is a workaround for cuda + gcc-7.2 compiler bug w/c++14
              const value_type *KOKKOS_RESTRICT ss = src + j * srcsize;
              /* */ value_type *KOKKOS_RESTRICT tt = tgt + j * srcsize;
              Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, srcsize),
                                   [&](const ordinal_type &i) { Kokkos::atomic_add(&tt[i], ss[i]); });
            });
        return;
      }
    }

    const ordinal_type *s_colidx = sbeg < send ? &info.gid_colidx(cur.gid_col_begin + srcbeg) : NULL;

    // loop over target
    // const size_type s2tsize = srcsize*sizeof(ordinal_type)*member.team_size();
    if constexpr(runOnHost) {
      ordinal_type *s2t = (ordinal_type *)buf;
      //const size_type s2tsize = srcsize * sizeof(ordinal_type);
      //TACHO_TEST_FOR_ABORT(bufsize < s2tsize, "bufsize is smaller than required s2t workspace");

      for (ordinal_type i = sbeg; i < send; ++i) {
        const ordinal_type tgtsid = info.sid_block_colidx(i).first;
        const auto &s = info.supernodes(tgtsid);
        {
          const ordinal_type tgtbeg = info.sid_block_colidx(s.sid_col_begin).second,
                             tgtend = info.sid_block_colidx(s.sid_col_end - 1).second, tgtsize = tgtend - tgtbeg;

          const ordinal_type *t_colidx = &info.gid_colidx(s.gid_col_begin + tgtbeg);
          for (ordinal_type k = 0, l = 0; k < srcsize; ++k) {
            s2t[k] = -1;
            for (; l < tgtsize && t_colidx[l] <= s_colidx[k]; ++l)
              if (s_colidx[k] == t_colidx[l]) {
                s2t[k] = l;
                break;
              }
          }
        }

        {
          UnmanagedViewType<value_type_matrix> U(s.u_buf, s.m, s.n);
          UnmanagedViewType<value_type_matrix> Lp(s.l_buf, s.n, s.m);
          const auto L = Kokkos::subview(Lp, range_type(s.m, s.n), Kokkos::ALL());

          ordinal_type ijbeg = 0;
          for (; s2t[ijbeg] == -1; ++ijbeg)
            ;

          // lock
          while (Kokkos::atomic_compare_exchange(&s.lock, 0, 1))
            TACHO_IMPL_PAUSE;
          Kokkos::store_fence();

          for (ordinal_type jj = ijbeg; jj < srcsize; ++jj) {
            const ordinal_type col = s2t[jj];
            for (ordinal_type ii = ijbeg; ii < srcsize; ++ii) {
              const ordinal_type row = s2t[ii];
              if (row < s.m) {
                U(row, col) += ABR(ii, jj);
                //Kokkos::atomic_add(&U(row, col), ABR(ii, jj));
                if (col >= s.m) {
                  L(col - s.m, row) += ABR(jj, ii);
                  //Kokkos::atomic_add(&L(col - s.m, row), ABR(jj, ii));
                }
              } else
                break;
            }
          }

          // unlock
          s.lock = 0;
          Kokkos::load_fence();
        }
      }
    } else {
      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(member, sbeg, send),
          [&, buf,
           srcsize](const ordinal_type &i) { // Value capture is a workaround for cuda + gcc-7.2 compiler bug w/c++14
            ordinal_type *s2t = ((ordinal_type *)(buf)) + member.team_rank() * srcsize;
            const auto &s = info.supernodes(info.sid_block_colidx(i).first);
            {
              const ordinal_type tgtbeg = info.sid_block_colidx(s.sid_col_begin).second,
                                 tgtend = info.sid_block_colidx(s.sid_col_end - 1).second, tgtsize = tgtend - tgtbeg;

              const ordinal_type *t_colidx = &info.gid_colidx(s.gid_col_begin + tgtbeg);
              Kokkos::parallel_for(
                  Kokkos::ThreadVectorRange(member, srcsize),
                  [&, t_colidx, s_colidx, tgtsize](
                      const ordinal_type &k) { // Value capture is a workaround for cuda + gcc-7.2 compiler bug w/c++14
                    s2t[k] = -1;
                    auto found = lower_bound(&t_colidx[0], &t_colidx[tgtsize - 1], s_colidx[k],
                                             [](ordinal_type left, ordinal_type right) { return left < right; });
                    if (s_colidx[k] == *found) {
                      s2t[k] = found - t_colidx;
                    }
                  });
            }
            {
              UnmanagedViewType<value_type_matrix> U(s.u_buf, s.m, s.n);
              UnmanagedViewType<value_type_matrix> Lp(s.l_buf, s.n, s.m);
              const auto L = Kokkos::subview(Lp, range_type(s.m, s.n), Kokkos::ALL());

              ordinal_type ijbeg = 0;
              for (; s2t[ijbeg] == -1; ++ijbeg)
                ;

              Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, ijbeg, srcsize), [&](const ordinal_type &ii) {
                const ordinal_type row = s2t[ii];
                if (row < s.m) {
                  for (ordinal_type jj = ijbeg; jj < srcsize; ++jj) {
                    const ordinal_type col = s2t[jj];
                    Kokkos::atomic_add(&U(row, col), ABR(ii, jj));
                    if (col >= s.m) {
                      Kokkos::atomic_add(&L(col - s.m, row), ABR(jj, ii));
                    }
                  }
                }
              });
            }
          });
    }
    return;
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void check(MemberType &member, supernode_type &s, const ordinal_type_array &fpiv) const {
    ordinal_type val(0);
    Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(member, s.m),
        [&](const ordinal_type &i, ordinal_type &update) {
          const ordinal_type fpiv_at_i = fpiv(i);
          update += (fpiv_at_i < 0 ? -fpiv_at_i : fpiv_at_i);
        },
        val);
    member.team_barrier();
    Kokkos::single(Kokkos::PerTeam(member), [&]() { s.do_not_apply_pivots = (val == 0); });
    return;
  }

  template <int Var> struct FactorizeTag {
    enum { variant = Var };
  };
  struct UpdateTag {};
  struct DummyTag {};

  template <typename MemberType, int Var>
  KOKKOS_INLINE_FUNCTION void operator()(const FactorizeTag<Var> &, const MemberType &member) const {
    const ordinal_type lid = member.league_rank();
    const ordinal_type p = _pbeg + lid;
    const ordinal_type sid = _level_sids(p);
    const ordinal_type mode = _compute_mode(sid);
    if (p < _pend && mode == 1) {
      using factorize_tag_type = FactorizeTag<Var>;

      const auto &s = _info.supernodes(sid);
      const ordinal_type m = s.m, n = s.n, n_m = n - m;
      const ordinal_type offm = s.row_begin;

      UnmanagedViewType<ordinal_type_array> P(_piv.data() + offm * 4, m * 4);

      const auto bufptr = _buf.data() + _buf_ptr(lid);
      if (factorize_tag_type::variant == 0) { // temporary to push to trilinos
        UnmanagedViewType<value_type_matrix> ABR(bufptr, n_m, n_m);
        factorize_var0(member, s, P, ABR);
      } else if (factorize_tag_type::variant == 1) {
        UnmanagedViewType<value_type_matrix> ABR(bufptr, n_m, n_m);
        UnmanagedViewType<value_type_matrix> T(bufptr, m, m);
        factorize_var1(member, s, P, T, ABR);
      } else if (factorize_tag_type::variant == 2 || factorize_tag_type::variant == 3) {
        UnmanagedViewType<value_type_matrix> ABR(bufptr, n_m, n_m);
        UnmanagedViewType<value_type_matrix> T(bufptr + ABR.span(), m, m);
        factorize_var2(member, s, P, T, ABR);
      }
    } else if (mode == -1) {
      Kokkos::printf("Error: TeamFunctorFactorizeChol, computing mode is not determined\n");
    } else {
      // skip
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const UpdateTag &, const MemberType &member) const {
    const ordinal_type lid = member.league_rank();
    const ordinal_type p = _pbeg + lid;
    if (p < _pend) {
      const ordinal_type sid = _level_sids(p);
      auto &s = _info.supernodes(sid);
      const ordinal_type n_m = s.n - s.m;
      value_type *bufptr = _buf.data() + _buf_ptr(lid);
      UnmanagedViewType<value_type_matrix> ABR(bufptr, n_m, n_m);
      update(member, s, ABR);

      const ordinal_type offm = s.row_begin;
      UnmanagedViewType<ordinal_type_array> fpiv(_piv.data() + offm * 4 + s.m, s.m);
      check(member, s, fpiv);
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
