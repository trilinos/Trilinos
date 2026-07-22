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
#ifndef __TACHO_CHOL_SUPERNODES_SERIAL_HPP__
#define __TACHO_CHOL_SUPERNODES_SERIAL_HPP__

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

template <> struct CholSupernodes<Algo::Workflow::Serial> {
  template <typename MemberType, typename SupernodeInfoType>
  KOKKOS_INLINE_FUNCTION static int factorize(MemberType &member, const SupernodeInfoType &info,
                                              const typename SupernodeInfoType::value_type_matrix &ABR,
                                              const ordinal_type sid) {
    using supernode_info_type = SupernodeInfoType;

    using value_type = typename supernode_info_type::value_type;
    using value_type_matrix = typename supernode_info_type::value_type_matrix;

    // algorithm choice
    using CholAlgoType = typename CholAlgorithm::type;
    using TrsmAlgoType = typename TrsmAlgorithm::type;
    using HerkAlgoType = typename HerkAlgorithm::type;

    // get current supernode
    const auto &s = info.supernodes(sid);

    // get panel pointer
    value_type *ptr = s.u_buf;

    // panel (s.m x s.n) is divided into ATL (m x m) and ATR (m x n)
    const ordinal_type m = s.m, n = s.n - s.m;

    // m and n are available, then factorize the supernode block
    if (m > 0) {
      UnmanagedViewType<value_type_matrix> ATL(ptr, m, m);
      ptr += m * m;
      Chol<Uplo::Upper, CholAlgoType>::invoke(member, ATL);

      if (n > 0) {
        const value_type one(1), zero(0);
        UnmanagedViewType<value_type_matrix> ATR(ptr, m, n); // ptr += m*n;
        Trsm<Side::Left, Uplo::Upper, Trans::ConjTranspose, TrsmAlgoType>::invoke(member, Diag::NonUnit(), one, ATL,
                                                                                  ATR);

        TACHO_TEST_FOR_ABORT(static_cast<ordinal_type>(ABR.extent(0)) != n ||
                                 static_cast<ordinal_type>(ABR.extent(1)) != n,
                             "ABR dimension does not match to supernodes");
        Herk<Uplo::Upper, Trans::ConjTranspose, HerkAlgoType>::invoke(member, -one, ATR, zero, ABR);
      }
    }
    return 0;
  }

  template <typename MemberType, typename SupernodeInfoType>
  KOKKOS_INLINE_FUNCTION static int update(MemberType &member, const SupernodeInfoType &info,
                                           const typename SupernodeInfoType::value_type_matrix &ABR,
                                           const ordinal_type sid, const size_type bufsize,
                                           /* */ void *buf, const bool update_lower = false) {
    using supernode_info_type = SupernodeInfoType;
    using value_type = typename supernode_info_type::value_type;
    using value_type_matrix = typename supernode_info_type::value_type_matrix;
    using range_type = typename supernode_info_type::range_type;

    const auto &cur = info.supernodes(sid);

    const ordinal_type sbeg = cur.sid_col_begin + 1, send = cur.sid_col_end - 1;

    const ordinal_type srcbeg = info.sid_block_colidx(sbeg).second, srcend = info.sid_block_colidx(send).second,
                       srcsize = srcend - srcbeg;

    static constexpr bool runOnHost = run_tacho_on_host_v <typename supernode_info_type::exec_space>;

    // short cut to direct update
    if ((send - sbeg) == 1) {
      const auto &s = info.supernodes(info.sid_block_colidx(sbeg).first);
      const ordinal_type tgtbeg = info.sid_block_colidx(s.sid_col_begin).second,
                         tgtend = info.sid_block_colidx(s.sid_col_end - 1).second, tgtsize = tgtend - tgtbeg;

      if (srcsize == tgtsize) {
        /* */ value_type *tgt = s.u_buf;
        const value_type *src = (value_type *)ABR.data();

        if constexpr(runOnHost) {
          // lock
          while (Kokkos::atomic_compare_exchange(&s.lock, 0, 1))
            TACHO_IMPL_PAUSE;
          Kokkos::store_fence();

          for (ordinal_type j = 0; j < srcsize; ++j) {
            const value_type *KOKKOS_RESTRICT ss = src + j * srcsize;
            /* */ value_type *KOKKOS_RESTRICT tt = tgt + j * srcsize;
            const ordinal_type iend = update_lower ? srcsize : j + 1;
            #if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
            #pragma unroll
            #endif
            for (ordinal_type i = 0; i < iend; ++i)
              tt[i] += ss[i];
          }

          // unlock
          s.lock = 0;
          Kokkos::load_fence();
        } else {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member, srcsize), [&](const ordinal_type &j) {
            const value_type *KOKKOS_RESTRICT ss = src + j * srcsize;
            /* */ value_type *KOKKOS_RESTRICT tt = tgt + j * srcsize;
            const ordinal_type iend = update_lower ? srcsize : j + 1;
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, iend),
                                 [&](const ordinal_type &i) { Kokkos::atomic_add(&tt[i], ss[i]); });
          });
        }

        return 0;
      }
    }

    const ordinal_type *s_colidx = sbeg < send ? &info.gid_colidx(cur.gid_col_begin + srcbeg) : NULL;

    // loop over target
    if constexpr(runOnHost) {
      ordinal_type *s2t = (ordinal_type *)buf;
      const size_type s2tsize = srcsize * sizeof(ordinal_type);
      TACHO_TEST_FOR_ABORT(bufsize < s2tsize, "bufsize is smaller than required s2t workspace");

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
            #if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
            #pragma unroll
            #endif
            for (ordinal_type ii = ijbeg; ii < srcsize; ++ii) {
              const ordinal_type row = s2t[ii];
              if (row < s.m) {
                U(row, col) += ABR(ii, jj);
                if (update_lower && col >= s.m) {
                  L(col - s.m, row) += ABR(jj, ii);
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
      // CUDA version
      const size_type s2tsize = srcsize * sizeof(ordinal_type) * member.team_size();
      TACHO_TEST_FOR_ABORT(bufsize < s2tsize, "bufsize is smaller than required s2t workspace");
#if 0 // single version for testing only 
        Kokkos::single(Kokkos::PerTeam(member), [&]() {
            ordinal_type *s2t = (ordinal_type*)buf;        
            for (ordinal_type i=sbeg;i<send;++i) {
              const auto &s = info.supernodes(info.sid_block_colidx(i).first);
              {
                const ordinal_type 
                  tgtbeg  = info.sid_block_colidx(s.sid_col_begin).second,
                  tgtend  = info.sid_block_colidx(s.sid_col_end-1).second,
                  tgtsize = tgtend - tgtbeg;
                
                const ordinal_type *t_colidx = &info.gid_colidx(s.gid_col_begin + tgtbeg);
                for (ordinal_type k=0,l=0;k<srcsize;++k) {
                  s2t[k] = -1;
                  for (;l<tgtsize && t_colidx[l] <= s_colidx[k];++l)
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

                ordinal_type ijbeg = 0; for (;s2t[ijbeg] == -1; ++ijbeg) ;
                
                for (ordinal_type jj=ijbeg;jj<srcsize;++jj) {
                  const ordinal_type col = s2t[jj];
                  for (ordinal_type ii=ijbeg;ii<srcsize;++ii) {
                    const ordinal_type row = s2t[ii];
                    if (row < s.m) {
                      Kokkos::atomic_add(&U(row, col), ABR(ii, jj));
                      if (update_lower && col >= s.m) {
                        Kokkos::atomic_add(&L(col-s.m, row), ABR(jj, ii));
                      }
                    } else break;
                  }
                }
              }
            }
          });
#else
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, sbeg, send), [&](const ordinal_type &i) {
        ordinal_type *s2t = ((ordinal_type *)(buf)) + member.team_rank() * srcsize;
        const auto &s = info.supernodes(info.sid_block_colidx(i).first);
        {
          const ordinal_type tgtbeg = info.sid_block_colidx(s.sid_col_begin).second,
                             tgtend = info.sid_block_colidx(s.sid_col_end - 1).second, tgtsize = tgtend - tgtbeg;

          const ordinal_type *t_colidx = &info.gid_colidx(s.gid_col_begin + tgtbeg);
          // for (ordinal_type k=0,l=0;k<srcsize;++k) {
          //   s2t[k] = -1;
          //   for (;l<tgtsize && t_colidx[l] <= s_colidx[k];++l)
          //     if (s_colidx[k] == t_colidx[l]) {
          //       s2t[k] = l;
          //       break;
          //     }
          // }
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, srcsize), [&](const ordinal_type &k) {
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
                if (update_lower && col >= s.m) {
                  Kokkos::atomic_add(&L(col - s.m, row), ABR(jj, ii));
                }
              }
            }
          });
          // Kokkos::parallel_for
          //   (Kokkos::ThreadVectorRange(member, srcsize-ijbeg), [&](const ordinal_type &jjj) {
          //     const ordinal_type jj = jjj + ijbeg;
          //     for (ordinal_type ii=ijbeg;ii<srcsize;++ii) {
          //       const ordinal_type row = s2t[ii];
          //       if (row < s.m)
          //         Kokkos::atomic_add(&A(row, s2t[jj]), ABR(ii, jj));
          //       else
          //         break;
          //     }
          //   });
          // const ordinal_type cnt = srcsize-ijbeg;
          // Kokkos::parallel_for
          //   (Kokkos::ThreadVectorRange(member, cnt*cnt), [&](const ordinal_type &idx) {
          //     const ordinal_type ii = idx%cnt + ijbeg;
          //     const ordinal_type jj = idx/cnt + ijbeg;
          //     const ordinal_type row = s2t[ii];
          //     const ordinal_type col = s2t[jj];
          //     if (row < s.m) Kokkos::atomic_add(&A(row, s2t[jj]), ABR(ii, jj));
          //   });
        }
      });
#endif
    }
    return 0;
  }

  template <typename MemberType, typename SupernodeInfoType>
  KOKKOS_INLINE_FUNCTION static int solve_lower(MemberType &member, const SupernodeInfoType &info,
                                                const typename SupernodeInfoType::value_type_matrix &xB,
                                                const ordinal_type sid) {
    using supernode_info_type = SupernodeInfoType;

    using value_type = typename supernode_info_type::value_type;
    using value_type_matrix = typename supernode_info_type::value_type_matrix;
    using range_type = typename supernode_info_type::range_type;

    const auto &s = info.supernodes(sid);

    using TrsmAlgoType = typename TrsmAlgorithm::type;
    using GemmAlgoType = typename GemmAlgorithm::type;

    // get panel pointer
    value_type *ptr = s.u_buf;

    // panel is divided into diagonal and interface block
    const ordinal_type m = s.m, n = s.n - s.m, nrhs = info.x.extent(1);

    // m and n are available, then factorize the supernode block
    if (m > 0) {
      const value_type one(1), zero(0);
      const ordinal_type offm = s.row_begin;
      UnmanagedViewType<value_type_matrix> AL(ptr, m, m);
      ptr += m * m;
      auto xT = Kokkos::subview(info.x, range_type(offm, offm + m), Kokkos::ALL());

      if (nrhs >= ThresholdSolvePhaseUsingBlas3)
        Trsm<Side::Left, Uplo::Upper, Trans::ConjTranspose, TrsmAlgoType>::invoke(member, Diag::NonUnit(), one, AL, xT);
      else
        Trsv<Uplo::Upper, Trans::ConjTranspose, TrsmAlgoType>::invoke(member, Diag::NonUnit(), AL, xT);

      if (n > 0) {
        UnmanagedViewType<value_type_matrix> AR(ptr, m, n); // ptr += m*n;
        if (nrhs >= ThresholdSolvePhaseUsingBlas3)
          Gemm<Trans::ConjTranspose, Trans::NoTranspose, GemmAlgoType>::invoke(member, -one, AR, xT, zero, xB);
        else
          Gemv<Trans::ConjTranspose, GemmAlgoType>::invoke(member, -one, AR, xT, zero, xB);
      }
    }
    return 0;
  }

  template <typename MemberType, typename SupernodeInfoType>
  KOKKOS_INLINE_FUNCTION static int update_solve_lower(MemberType &member, const SupernodeInfoType &info,
                                                       const typename SupernodeInfoType::value_type_matrix &xB,
                                                       const ordinal_type sid) {

    static constexpr bool runOnHost = run_tacho_on_host_v <typename SupernodeInfoType::exec_space>;

    const auto &cur = info.supernodes(sid);
    const ordinal_type sbeg = cur.sid_col_begin + 1, send = cur.sid_col_end - 1;

    const ordinal_type m = xB.extent(0), n = xB.extent(1);
    TACHO_TEST_FOR_ABORT(m != (cur.n - cur.m), "# of rows in xB does not match to super blocksize in sid");

    if constexpr(runOnHost) {
      for (ordinal_type i = sbeg, is = 0; i < send; ++i) {
        const ordinal_type tbeg = info.sid_block_colidx(i).second, tend = info.sid_block_colidx(i + 1).second;

        // lock
        const auto &s = info.supernodes(info.sid_block_colidx(i).first);
        while (Kokkos::atomic_compare_exchange(&s.lock, 0, 1))
          TACHO_IMPL_PAUSE;
        Kokkos::store_fence();

        // both src and tgt increase index
        for (ordinal_type it = tbeg; it < tend; ++it, ++is) {
          const ordinal_type row = info.gid_colidx(cur.gid_col_begin + it);
          for (ordinal_type j = 0; j < n; ++j)
            info.x(row, j) += xB(is, j);
        }

        // unlock
        s.lock = 0;
        Kokkos::load_fence();
      }
    } else {
      Kokkos::single(Kokkos::PerTeam(member), [&]() {
        for (ordinal_type i = sbeg, is = 0; i < send; ++i) {
          const ordinal_type tbeg = info.sid_block_colidx(i).second, tend = info.sid_block_colidx(i + 1).second;

          for (ordinal_type it = tbeg; it < tend; ++it, ++is) {
            const ordinal_type row = info.gid_colidx(cur.gid_col_begin + it);
            for (ordinal_type j = 0; j < n; ++j)
              Kokkos::atomic_add(&info.x(row, j), xB(is, j));
          }
        }
      });
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

    using GemmAlgoType = typename GemmAlgorithm::type;
    using TrsmAlgoType = typename TrsmAlgorithm::type;

    // get current supernode
    const auto &s = info.supernodes(sid);

    // get supernode panel pointer
    value_type *ptr = s.u_buf;

    // panel is divided into diagonal and interface block
    const ordinal_type m = s.m, n = s.n - s.m, nrhs = info.x.extent(1);

    // m and n are available, then factorize the supernode block
    if (m > 0) {
      const value_type one(1);
      const UnmanagedViewType<value_type_matrix> AL(ptr, m, m);
      ptr += m * m;

      const ordinal_type offm = s.row_begin;
      const auto xT = Kokkos::subview(info.x, range_type(offm, offm + m), Kokkos::ALL());

      if (n > 0) {
        const UnmanagedViewType<value_type_matrix> AR(ptr, m, n); // ptr += m*n;
        if (nrhs >= ThresholdSolvePhaseUsingBlas3)
          Gemm<Trans::NoTranspose, Trans::NoTranspose, GemmAlgoType>::invoke(member, -one, AR, xB, one, xT);
        else
          Gemv<Trans::NoTranspose, GemmAlgoType>::invoke(member, -one, AR, xB, one, xT);
      }
      if (nrhs >= ThresholdSolvePhaseUsingBlas3)
        Trsm<Side::Left, Uplo::Upper, Trans::NoTranspose, TrsmAlgoType>::invoke(member, Diag::NonUnit(), one, AL, xT);
      else
        Trsv<Uplo::Upper, Trans::NoTranspose, TrsmAlgoType>::invoke(member, Diag::NonUnit(), AL, xT);
    }
    return 0;
  }

  template <typename MemberType, typename SupernodeInfoType>
  KOKKOS_INLINE_FUNCTION static int update_solve_upper(MemberType &member, const SupernodeInfoType &info,
                                                       const typename SupernodeInfoType::value_type_matrix &xB,
                                                       const ordinal_type sid) {

    static constexpr bool runOnHost = run_tacho_on_host_v <typename SupernodeInfoType::exec_space>;

    const auto &s = info.supernodes(sid);

    const ordinal_type m = xB.extent(0), n = xB.extent(1);
    TACHO_TEST_FOR_ABORT(m != (s.n - s.m), "# of rows in xB does not match to super blocksize in sid");

    const ordinal_type goffset = s.gid_col_begin + s.m;
    if constexpr(runOnHost) {
      for (ordinal_type j = 0; j < n; ++j)
        for (ordinal_type i = 0; i < m; ++i) {
          const ordinal_type row = info.gid_colidx(i + goffset);
          xB(i, j) = info.x(row, j);
        }
    } else {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n), [&](const ordinal_type &j) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, m), [&](const ordinal_type &i) {
          const ordinal_type row = info.gid_colidx(i + goffset);
          xB(i, j) = info.x(row, j);
        });
      });
    }

    return 0;
  }

  template <typename MemberType, typename SupernodeInfoType>
  KOKKOS_INLINE_FUNCTION static int
  factorize_recursive_serial(MemberType &member, const SupernodeInfoType &info, const ordinal_type sid,
                             const bool final, typename SupernodeInfoType::value_type *buf, const size_type bufsize) {
    using supernode_info_type = SupernodeInfoType;

    using value_type = typename supernode_info_type::value_type;
    using value_type_matrix = typename supernode_info_type::value_type_matrix;

    static constexpr bool runOnHost = run_tacho_on_host_v <typename supernode_info_type::exec_space>;

    const auto &s = info.supernodes(sid);

    if (final) {
      // serial recursion
      for (ordinal_type i = 0; i < s.nchildren; ++i)
        factorize_recursive_serial(member, info, s.children[i], final, buf, bufsize);
    }

    {
      const ordinal_type n = s.n - s.m;

      size_type bufsize_required;
      if constexpr(runOnHost)
        bufsize_required = n * (n + 1) * sizeof(value_type);
      else
        bufsize_required = n * (n + member.team_size()) * sizeof(value_type);

      TACHO_TEST_FOR_ABORT(bufsize < bufsize_required, "bufsize is smaller than required");

      UnmanagedViewType<value_type_matrix> ABR((value_type *)buf, n, n);

      CholSupernodes<Algo::Workflow::Serial>::factorize(member, info, ABR, sid);

      CholSupernodes<Algo::Workflow::Serial>::update(member, info, ABR, sid, bufsize - ABR.span() * sizeof(value_type),
                                                     (void *)((value_type *)buf + ABR.span()));
    }
    return 0;
  }

  template <typename MemberType, typename SupernodeInfoType>
  KOKKOS_INLINE_FUNCTION static int
  solve_lower_recursive_serial(MemberType &member, const SupernodeInfoType &info, const ordinal_type sid,
                               const bool final, typename SupernodeInfoType::value_type *buf, const size_type bufsize) {
    using supernode_info_type = SupernodeInfoType;

    using value_type = typename supernode_info_type::value_type;
    using value_type_matrix = typename supernode_info_type::value_type_matrix;

    const auto &s = info.supernodes(sid);

    if (final) {
      // serial recursion
      for (ordinal_type i = 0; i < s.nchildren; ++i)
        solve_lower_recursive_serial(member, info, s.children[i], final, buf, bufsize);
    }

    {
      const ordinal_type n = s.n - s.m;
      const ordinal_type nrhs = info.x.extent(1);
      const size_type bufsize_required = n * nrhs * sizeof(value_type);

      TACHO_TEST_FOR_ABORT(bufsize < bufsize_required, "bufsize is smaller than required");

      UnmanagedViewType<value_type_matrix> xB((value_type *)buf, n, nrhs);

      CholSupernodes<Algo::Workflow::Serial>::solve_lower(member, info, xB, sid);

      CholSupernodes<Algo::Workflow::Serial>::update_solve_lower(member, info, xB, sid);
    }
    return 0;
  }

  template <typename MemberType, typename SupernodeInfoType>
  KOKKOS_INLINE_FUNCTION static int solve_upper_recursive_serial(MemberType &member, const SupernodeInfoType &info,
                                                                 const ordinal_type sid, const bool final,
                                                                 typename SupernodeInfoType::value_type *buf,
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

      CholSupernodes<Algo::Workflow::Serial>::solve_upper(member, info, xB, sid);
    }

    if (final) {
      // serial recursion
      for (ordinal_type i = 0; i < s.nchildren; ++i)
        solve_upper_recursive_serial(member, info, s.children[i], final, buf, bufsize);
    }
    return 0;
  }
};
} // namespace Tacho

#endif
