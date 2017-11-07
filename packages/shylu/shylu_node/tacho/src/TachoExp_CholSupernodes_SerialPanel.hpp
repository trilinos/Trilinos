#ifndef __TACHOEXP_CHOL_SUPERNODES_SERIAL_PANEL_HPP__
#define __TACHOEXP_CHOL_SUPERNODES_SERIAL_PANEL_HPP__

/// \file TachoExp_CholSupernodes.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

namespace Tacho {

  namespace Experimental {

    template<>
    struct CholSupernodes<Algo::Workflow::SerialPanel> {
      template<typename SchedulerType,
               typename MemberType,
               typename SupernodeInfoType>
      KOKKOS_INLINE_FUNCTION
      static int
      factorize(const SchedulerType &sched,
                const MemberType &member,
                const SupernodeInfoType &info,
                const ordinal_type sid) {
        typedef SupernodeInfoType supernode_info_type;

        typedef typename supernode_info_type::value_type value_type;
        typedef typename supernode_info_type::value_type_matrix value_type_matrix;

        // get current supernode 
        const auto &s = info.supernodes(sid);

        // get panel pointer
        value_type *ptr = s.buf;

        // panel (s.m x s.n) is divided into ATL (m x m) and ATR (m x n)
        const ordinal_type m = s.m, n = s.n - s.m;

        // m is available, then factorize the supernode block
        if (m > 0) {
          UnmanagedViewType<value_type_matrix> ATL(ptr, m, m); ptr += m*m;
          Chol<Uplo::Upper,Algo::External>::invoke(sched, member, ATL);

          // n is available, then solve interface block
          if (n > 0) {
            UnmanagedViewType<value_type_matrix> ATR(ptr, m, n);
            Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::External>
              ::invoke(sched, member, Diag::NonUnit(), 1.0, ATL, ATR);
          }
        }
        return 0;
      }

      template<typename SchedulerType,
               typename MemberType,
               typename SupernodeInfoType>
      KOKKOS_INLINE_FUNCTION
      static int
      update(const SchedulerType &sched,
             const MemberType &member,
             const SupernodeInfoType &info,
             const ordinal_type offn, // ATR and ABR panel offset
             const ordinal_type np,   // ATR and ABR panel width
             const ordinal_type sid,
             const size_type bufsize, // ABR size + additional
             /* */ void *buf,
             const bool use_atomic = true) {
        typedef SupernodeInfoType supernode_info_type;

        typedef typename supernode_info_type::value_type value_type;
        typedef typename supernode_info_type::value_type_matrix value_type_matrix;
        typedef typename supernode_info_type::dense_block_type dense_block_type;

        // get current supernode
        const auto &cur = info.supernodes(sid);

        // panel (cur.m x cur.n) is divided into ATL (m x m) and ATR (m x n)
        const ordinal_type 
          m = cur.m, n = cur.n - cur.m , 
          nb = min(np, n - offn), nn = offn + nb;

        // m and n are available, then factorize the supernode block
        if (m > 0 && n > 0) {
          // ** update
          const ordinal_type 
            sbeg = cur.sid_col_begin + 1, send = cur.sid_col_end - 1;
          
          const ordinal_type 
            srcbeg  = info.sid_block_colidx(sbeg).second, 
            srcend  = info.sid_block_colidx(send).second, 
            srcsize = srcend - srcbeg;

          TACHO_TEST_FOR_ABORT(bufsize < (srcsize*sizeof(ordinal_type) + nn*nb*sizeof(value_type)), 
                               "bufsize is smaller than required workspace");        
          
          UnmanagedViewType<value_type_matrix> ABL(cur.buf    +    m*m, m, nn);
          UnmanagedViewType<value_type_matrix> ATR(ABL.data() + offn*m, m, nb);

          value_type *ptr = (value_type*)buf;          
          UnmanagedViewType<value_type_matrix> ABR(ptr, nn, nb); ptr += ABR.span();

          if (offn == 0 && nb == n)
            Herk<Uplo::Upper,Trans::ConjTranspose,Algo::External>
              ::invoke(sched, member, -1.0, ATR, 0.0, ABR);            
          else 
            Gemm<Trans::ConjTranspose,Trans::NoTranspose,Algo::External>
              ::invoke(sched, member, -1.0, ABL, ATR, 0.0, ABR);  
          
          // short cut to direct update
          if ((send - sbeg) == 1) {
            const auto &s = info.supernodes(info.sid_block_colidx(sbeg).first);
            const ordinal_type 
              tgtbeg  = info.sid_block_colidx(s.sid_col_begin).second,
              tgtend  = info.sid_block_colidx(s.sid_col_end-1).second,
              tgtsize = tgtend - tgtbeg;
            
            if (srcsize == tgtsize) {
              /* */ value_type *tgt = s.buf;
              const value_type *src = (value_type*)ABR.data();
              
              // lock
              while (Kokkos::atomic_compare_exchange(&s.lock, 0, 1)) KOKKOS_IMPL_PAUSE;
              Kokkos::store_fence();            
              
              for (ordinal_type js=0;js<nb;++js) {
                const ordinal_type jt = js + offn;
                const value_type *__restrict__ ss = src + js*srcsize;
                /* */ value_type *__restrict__ tt = tgt + jt*srcsize;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
                for (ordinal_type i=0;i<=jt;++i)
                  tt[i] += ss[i];
              }
              
              // unlock
              s.lock = 0;
              Kokkos::load_fence();
              
              return 0;
            }
          } 
          
          dense_block_type A;
          ordinal_type *s2t = (ordinal_type*)ptr;
          
          // loop over target
          const ordinal_type *s_colidx = sbeg < send ? &info.gid_colidx(cur.gid_col_begin + srcbeg) : NULL;
          for (ordinal_type i=sbeg;i<send;++i) {
            const auto &s = info.supernodes(info.sid_block_colidx(i).first);
            {
              const ordinal_type 
                tgtbeg  = info.sid_block_colidx(s.sid_col_begin).second,
                tgtend  = info.sid_block_colidx(s.sid_col_end-1).second,
                tgtsize = tgtend - tgtbeg;
              
              const ordinal_type *t_colidx = &info.gid_colidx(s.gid_col_begin + tgtbeg);
              for (ordinal_type k=0,l=0;k<nn;++k) {
                s2t[k] = -1;
                for (;l<tgtsize && t_colidx[l] <= s_colidx[k];++l)
                  if (s_colidx[k] == t_colidx[l]) {
                    s2t[k] = l; 
                    break;
                  }
              }
            }
            
            {
              A.set_view(s.m, s.n);
              A.attach_buffer(1, s.m, s.buf);
              
              ordinal_type ijbeg = 0; for (;s2t[ijbeg] == -1; ++ijbeg) ;
              
              // lock
              if (use_atomic) {
                for (ordinal_type jj=max(ijbeg,offn);jj<nn;++jj) {
                  const ordinal_type js = jj - offn;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
                  for (ordinal_type ii=ijbeg;ii<nn;++ii) {
                    const ordinal_type row = s2t[ii];
                    if (row < s.m) Kokkos::atomic_fetch_add(&A(row, s2t[jj]), ABR(ii, js));
                    else break;
                  }
                }
              } else {
                while (Kokkos::atomic_compare_exchange(&s.lock, 0, 1)) KOKKOS_IMPL_PAUSE;
                Kokkos::store_fence();            
                
                for (ordinal_type jj=max(ijbeg,offn);jj<nn;++jj) {
                  const ordinal_type js = jj - offn;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
                  for (ordinal_type ii=ijbeg;ii<nn;++ii) {
                    const ordinal_type row = s2t[ii];
                    if (row < s.m) A(row, s2t[jj]) += ABR(ii, js);
                    else break;
                  }
                }
                // unlock
                s.lock = 0;
                Kokkos::load_fence();
              }
            }
          }
        }
        return 0;
      }

      template<typename SchedulerType,
               typename MemberType,
               typename SupernodeInfoType>
      KOKKOS_INLINE_FUNCTION
      static int
      factorize_recursive_serial(const SchedulerType &sched,
                                 const MemberType &member,
                                 const SupernodeInfoType &info,
                                 const ordinal_type sid,
                                 const bool final,
                                 typename SupernodeInfoType::value_type *buf,
                                 const ordinal_type bufsize,
                                 const ordinal_type np) {
        typedef SupernodeInfoType supernode_info_type;

        typedef typename supernode_info_type::value_type value_type;
        //typedef typename supernode_info_type::value_type_matrix value_type_matrix;

        const auto &s = info.supernodes(sid);

        if (final) {
          // serial recursion
          for (ordinal_type i=0;i<s.nchildren;++i)
            factorize_recursive_serial(sched, member, info, 
                                       s.children[i], final, buf, bufsize,
                                       np);
        }

        {
          const size_type n = s.n - s.m, bufsize_required = n*(min(np,n)+1)*sizeof(value_type);

          TACHO_TEST_FOR_ABORT(bufsize < static_cast<ordinal_type>(bufsize_required), 
                               "bufsize is smaller than required");

          CholSupernodes<Algo::Workflow::SerialPanel>
            ::factorize(sched, member, info, sid);

          for (ordinal_type offn=0; offn<static_cast<ordinal_type>(n); offn+=np) {
            CholSupernodes<Algo::Workflow::SerialPanel>
              ::update(sched, member, info, offn, np, sid, bufsize, (void*)buf);
          }
        }
        return 0;
      }

    };
  }
}

#endif
