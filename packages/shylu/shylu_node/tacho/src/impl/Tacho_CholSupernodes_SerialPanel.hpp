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
#ifndef __TACHO_CHOL_SUPERNODES_SERIAL_PANEL_HPP__
#define __TACHO_CHOL_SUPERNODES_SERIAL_PANEL_HPP__

/// \file Tacho_CholSupernodes.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

#include "Tacho_Lapack_External.hpp"
#include "Tacho_Lapack_Team.hpp"

#include "Tacho_Blas_External.hpp"
#include "Tacho_Blas_Team.hpp"

#include "Tacho_Chol.hpp"
#include "Tacho_Chol_External.hpp"
#include "Tacho_Chol_Internal.hpp"

#include "Tacho_Trsm.hpp"
#include "Tacho_Trsm_External.hpp"
#include "Tacho_Trsm_Internal.hpp"

#include "Tacho_Herk.hpp"
#include "Tacho_Herk_External.hpp"
#include "Tacho_Herk_Internal.hpp"

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

template <> struct CholSupernodes<Algo::Workflow::SerialPanel> {
  template <typename MemberType, typename SupernodeInfoType>
  KOKKOS_INLINE_FUNCTION static int factorize(MemberType &member, const SupernodeInfoType &info,
                                              const ordinal_type sid) {
    typedef SupernodeInfoType supernode_info_type;

    typedef typename supernode_info_type::value_type value_type;
    typedef typename supernode_info_type::value_type_matrix value_type_matrix;

    // algorithm choice
    using CholAlgoType = typename CholAlgorithm::type;
    using TrsmAlgoType = typename TrsmAlgorithm::type;

    // get current supernode
    const auto &s = info.supernodes(sid);

    // get panel pointer
    value_type *ptr = s.buf;

    // panel (s.m x s.n) is divided into ATL (m x m) and ATR (m x n)
    const ordinal_type m = s.m, n = s.n - s.m;

    // m is available, then factorize the supernode block
    if (m > 0) {
      UnmanagedViewType<value_type_matrix> ATL(ptr, m, m);
      ptr += m * m;
      Chol<Uplo::Upper, CholAlgoType>::invoke(member, ATL);

      // n is available, then solve interface block
      if (n > 0) {
        UnmanagedViewType<value_type_matrix> ATR(ptr, m, n);
        Trsm<Side::Left, Uplo::Upper, Trans::ConjTranspose, TrsmAlgoType>::invoke(member, Diag::NonUnit(), 1.0, ATL,
                                                                                  ATR);
      }
    }
    return 0;
  }

  template <typename MemberType, typename SupernodeInfoType>
  KOKKOS_INLINE_FUNCTION static int update(MemberType &member, const SupernodeInfoType &info,
                                           const ordinal_type offn, // ATR and ABR panel offset
                                           const ordinal_type np,   // ATR and ABR panel width
                                           const ordinal_type sid,
                                           const size_type bufsize, // ABR size + additional
                                           /* */ void *buf) {
    typedef SupernodeInfoType supernode_info_type;

    typedef typename supernode_info_type::value_type value_type;
    typedef typename supernode_info_type::value_type_matrix value_type_matrix;

    static constexpr bool runOnHost = run_tacho_on_host_v<typename supernode_info_type::exec_space>;

    // algorithm choice
    using HerkAlgoType = typename HerkAlgorithm::type;
    using GemmAlgoType = typename GemmAlgorithm::type;

    if constexpr(!runOnHost) {
      member.team_barrier();
    }
    // get current supernode
    const auto &cur = info.supernodes(sid);

    // panel (cur.m x cur.n) is divided into ATL (m x m) and ATR (m x n)
    const ordinal_type m = cur.m, n = cur.n - cur.m, nb = min(np, n - offn), nn = offn + nb;

    // m and n are available, then factorize the supernode block
    if (m > 0 && n > 0) {
      // ** update
      const ordinal_type sbeg = cur.sid_col_begin + 1, send = cur.sid_col_end - 1;

      const ordinal_type srcbeg = info.sid_block_colidx(sbeg).second, srcend = info.sid_block_colidx(send).second,
                         srcsize = srcend - srcbeg;

      if constexpr(runOnHost) {
        TACHO_TEST_FOR_ABORT(bufsize < size_type(srcsize * sizeof(ordinal_type) + nn * nb * sizeof(value_type)),
                             "bufsize is smaller than required workspace");
      } else {
        TACHO_TEST_FOR_ABORT(
            bufsize < size_type(srcsize * sizeof(ordinal_type) * member.team_size() + nn * nb * sizeof(value_type)),
            "bufsize is smaller than required workspace");
      }

      UnmanagedViewType<value_type_matrix> ABL(cur.buf + m * m, m, nn);
      UnmanagedViewType<value_type_matrix> ATR(ABL.data() + offn * m, m, nb);

      value_type *ptr = (value_type *)buf;
      UnmanagedViewType<value_type_matrix> ABR(ptr, nn, nb);
      ptr += ABR.span();

      if (offn == 0 && nb == n)
        Herk<Uplo::Upper, Trans::ConjTranspose, HerkAlgoType>::invoke(member, -1.0, ATR, 0.0, ABR);
      else
        Gemm<Trans::ConjTranspose, Trans::NoTranspose, GemmAlgoType>::invoke(member, -1.0, ABL, ATR, 0.0, ABR);

      // short cut to direct update
      if ((send - sbeg) == 1) {
        const auto &s = info.supernodes(info.sid_block_colidx(sbeg).first);
        const ordinal_type tgtbeg = info.sid_block_colidx(s.sid_col_begin).second,
                           tgtend = info.sid_block_colidx(s.sid_col_end - 1).second, tgtsize = tgtend - tgtbeg;

        if (srcsize == tgtsize) {
          /* */ value_type *tgt = s.buf;
          const value_type *src = (value_type *)ABR.data();

          if constexpr(runOnHost) {
            switch (info.front_update_mode) {
            case 1: {
              for (ordinal_type js = 0; js < nb; ++js) {
                const ordinal_type jt = js + offn;
                const value_type *KOKKOS_RESTRICT ss = src + js * srcsize;
                /* */ value_type *KOKKOS_RESTRICT tt = tgt + jt * srcsize;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
                for (ordinal_type i = 0; i <= jt; ++i)
                  Kokkos::atomic_fetch_add(&tt[i], ss[i]);
              }
              break;
            }
            case 0: {
              // lock
              while (Kokkos::atomic_compare_exchange(&s.lock, 0, 1))
                TACHO_IMPL_PAUSE;
              Kokkos::store_fence();

              for (ordinal_type js = 0; js < nb; ++js) {
                const ordinal_type jt = js + offn;
                const value_type *KOKKOS_RESTRICT ss = src + js * srcsize;
                /* */ value_type *KOKKOS_RESTRICT tt = tgt + jt * srcsize;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
                for (ordinal_type i = 0; i <= jt; ++i)
                  tt[i] += ss[i];
              }

              // unlock
              s.lock = 0;
              Kokkos::load_fence();
              break;
            }
            }
          } else {
            Kokkos::parallel_for(Kokkos::TeamThreadRange(member, nb), [&](const ordinal_type &js) {
              const ordinal_type jt = js + offn;
              const value_type *KOKKOS_RESTRICT ss = src + js * srcsize;
              /* */ value_type *KOKKOS_RESTRICT tt = tgt + jt * srcsize;
              Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, jt + 1),
                                   [&](const ordinal_type &i) { Kokkos::atomic_fetch_add(&tt[i], ss[i]); });
            });
          }

          return 0;
        }
      }

      if constexpr(runOnHost) {
        // loop over target
        const ordinal_type *s_colidx = sbeg < send ? &info.gid_colidx(cur.gid_col_begin + srcbeg) : NULL;
        for (ordinal_type i = sbeg; i < send; ++i) {
          ordinal_type *s2t = (ordinal_type *)ptr;
          const auto &s = info.supernodes(info.sid_block_colidx(i).first);
          {
            const ordinal_type tgtbeg = info.sid_block_colidx(s.sid_col_begin).second,
                               tgtend = info.sid_block_colidx(s.sid_col_end - 1).second, tgtsize = tgtend - tgtbeg;

            const ordinal_type *t_colidx = &info.gid_colidx(s.gid_col_begin + tgtbeg);
            for (ordinal_type k = 0, l = 0; k < nn; ++k) {
              s2t[k] = -1;
              for (; l < tgtsize && t_colidx[l] <= s_colidx[k]; ++l)
                if (s_colidx[k] == t_colidx[l]) {
                  s2t[k] = l;
                  break;
                }
            }
          }

          {
            UnmanagedViewType<value_type_matrix> A(s.u_buf, s.m, s.n);

            ordinal_type ijbeg = 0;
            for (; s2t[ijbeg] == -1; ++ijbeg)
              ;

            // lock
            switch (info.front_update_mode) {
            case 1: {
              for (ordinal_type jj = max(ijbeg, offn); jj < nn; ++jj) {
                const ordinal_type js = jj - offn;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
                for (ordinal_type ii = ijbeg; ii < nn; ++ii) {
                  const ordinal_type row = s2t[ii];
                  if (row < s.m)
                    Kokkos::atomic_fetch_add(&A(row, s2t[jj]), ABR(ii, js));
                  else
                    break;
                }
              }
              break;
            }
            case 0: {
              while (Kokkos::atomic_compare_exchange(&s.lock, 0, 1))
                TACHO_IMPL_PAUSE;
              Kokkos::store_fence();

              for (ordinal_type jj = max(ijbeg, offn); jj < nn; ++jj) {
                const ordinal_type js = jj - offn;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
                for (ordinal_type ii = ijbeg; ii < nn; ++ii) {
                  const ordinal_type row = s2t[ii];
                  if (row < s.m)
                    A(row, s2t[jj]) += ABR(ii, js);
                  else
                    break;
                }
              }
              // unlock
              s.lock = 0;
              Kokkos::load_fence();
              break;
            }
            }
          }
        }
      } else {
#if 1
        const ordinal_type *s_colidx = sbeg < send ? &info.gid_colidx(cur.gid_col_begin + srcbeg) : NULL;
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, sbeg, send), [&](const ordinal_type &i) {
          ordinal_type *s2t = (((ordinal_type *)ptr) + (member.team_rank() * nn));
          const auto &s = info.supernodes(info.sid_block_colidx(i).first);
          {
            const ordinal_type tgtbeg = info.sid_block_colidx(s.sid_col_begin).second,
                               tgtend = info.sid_block_colidx(s.sid_col_end - 1).second, tgtsize = tgtend - tgtbeg;

            const ordinal_type *t_colidx = &info.gid_colidx(s.gid_col_begin + tgtbeg);
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, nn), [&](const ordinal_type &k) {
              s2t[k] = -1;
              auto found = lower_bound(&t_colidx[0], &t_colidx[tgtsize - 1], s_colidx[k],
                                       [](ordinal_type left, ordinal_type right) { return left < right; });
              if (s_colidx[k] == *found) {
                s2t[k] = found - t_colidx;
              }
            });
          }

          {
            UnmanagedViewType<value_type_matrix> A(s.u_buf, s.m, s.n);

            ordinal_type ijbeg = 0;
            for (; s2t[ijbeg] == -1; ++ijbeg)
              ;

            Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, nn - ijbeg), [&](const ordinal_type &iii) {
              const ordinal_type ii = ijbeg + iii;
              const ordinal_type row = s2t[ii];
              if (row < s.m)
                for (ordinal_type jj = max(ijbeg, offn); jj < nn; ++jj) {
                  const ordinal_type js = jj - offn;
                  Kokkos::atomic_add(&A(row, s2t[jj]), ABR(ii, js));
                }
            });
            // for (ordinal_type jj=max(ijbeg,offn);jj<nn;++jj) {
            //   const ordinal_type js = jj - offn;
            //   Kokkos::parallel_for
            //     (Kokkos::ThreadVectorRange(member, nn - ijbeg), [&](const ordinal_type &iii) {
            //       const ordinal_type ii = iii + ijbeg;
            //       const ordinal_type row = s2t[ii];
            //       if (row < s.m) Kokkos::atomic_fetch_add(&A(row, s2t[jj]), ABR(ii, js));
            //     });
            // }
            // const ordinal_type ijtmp = max(ijbeg,offn);
            // Kokkos::parallel_for
            //   (Kokkos::ThreadVectorRange(member, nn - ijtmp), [&](const ordinal_type &jjj) {
            //     const ordinal_type jj = jjj + ijtmp;
            //     const ordinal_type js = jj - offn;
            //     for (ordinal_type ii=ijbeg;ii<nn;++ii) {
            //       const ordinal_type row = s2t[ii];
            //       if (row < s.m)
            //         Kokkos::atomic_fetch_add(&A(row, s2t[jj]), ABR(ii, js));
            //       else
            //         break;
            //     }
            //   });
          }
        });
#endif
      }
    }
    return 0;
  }

  template <typename MemberType, typename SupernodeInfoType>
  KOKKOS_INLINE_FUNCTION static int factorize_recursive_serial(MemberType &member, const SupernodeInfoType &info,
                                                               const ordinal_type sid, const bool final,
                                                               typename SupernodeInfoType::value_type *buf,
                                                               const size_type bufsize, const size_type np) {
    typedef SupernodeInfoType supernode_info_type;

    typedef typename supernode_info_type::value_type value_type;

    static constexpr bool runOnHost = run_tacho_on_host_v<typename supernode_info_type::exec_space>;

    const auto &s = info.supernodes(sid);

    if (final) {
      // serial recursion
      for (ordinal_type i = 0; i < s.nchildren; ++i)
        factorize_recursive_serial(member, info, s.children[i], final, buf, bufsize, np);
    }

    {
      const ordinal_type n = s.n - s.m;

      size_type bufsize_required;
      if constexpr(runOnHost)
        bufsize_required = n * (min(np, n) + 1) * sizeof(value_type);
      else
        bufsize_required = n * (min(np, n) + member.team_size()) * sizeof(value_type);

      TACHO_TEST_FOR_ABORT(bufsize < bufsize_required, "bufsize is smaller than required");

      CholSupernodes<Algo::Workflow::SerialPanel>::factorize(member, info, sid);

      for (ordinal_type offn = 0; offn < n; offn += np) {
        CholSupernodes<Algo::Workflow::SerialPanel>::update(member, info, offn, np, sid, bufsize, (void *)buf);
      }
    }
    return 0;
  }
};
} // namespace Tacho

#endif
