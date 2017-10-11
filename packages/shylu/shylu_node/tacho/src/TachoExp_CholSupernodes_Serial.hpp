#ifndef __TACHOEXP_CHOL_SUPERNODES_SERIAL_HPP__
#define __TACHOEXP_CHOL_SUPERNODES_SERIAL_HPP__

/// \file TachoExp_CholSupernodes.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

namespace Tacho {

  namespace Experimental {

    template<>
    struct CholSupernodes<Algo::Workflow::Serial> {
      template<typename SchedulerType,
               typename MemberType,
               typename SupernodeInfoType>
      KOKKOS_INLINE_FUNCTION
      static int
      factorize(const SchedulerType &sched,
                const MemberType &member,
                const SupernodeInfoType &info,
                const typename SupernodeInfoType::value_type_matrix &ABR,
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

        // m and n are available, then factorize the supernode block
        if (m > 0) {
          UnmanagedViewType<value_type_matrix> ATL(ptr, m, m); ptr += m*m;
          Chol<Uplo::Upper,Algo::External>
            ::invoke(sched, member, ATL);

          if (n > 0) {
            UnmanagedViewType<value_type_matrix> ATR(ptr, m, n); // ptr += m*n;
            Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::External>
              ::invoke(sched, member, Diag::NonUnit(), 1.0, ATL, ATR);

            TACHO_TEST_FOR_ABORT(static_cast<ordinal_type>(ABR.dimension_0()) != n ||
                                 static_cast<ordinal_type>(ABR.dimension_1()) != n,
                                 "ABR dimension does not match to supernodes");
            Herk<Uplo::Upper,Trans::ConjTranspose,Algo::External>
              ::invoke(sched, member, -1.0, ATR, 0.0, ABR);
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
             const typename SupernodeInfoType::value_type_matrix &ABR,
             const ordinal_type sid,
             const size_type bufsize,
             /* */ void *buf) {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
#if defined( TACHO_PROFILE_TIME_PER_THREAD )
        Kokkos::Impl::Timer timer;
        timer.reset();
#endif
#endif
        typedef SupernodeInfoType supernode_info_type;

        typedef typename supernode_info_type::value_type value_type;
        //typedef typename supernode_info_type::ordinal_type_array ordinal_type_array;
        typedef typename supernode_info_type::dense_block_type dense_block_type;

        const auto &cur = info.supernodes(sid);

        const ordinal_type 
          sbeg = cur.sid_col_begin + 1, send = cur.sid_col_end - 1;

        const ordinal_type 
          srcbeg  = info.sid_block_colidx(sbeg).second, 
          srcend  = info.sid_block_colidx(send).second, 
          srcsize = srcend - srcbeg;
        
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

            for (ordinal_type j=0;j<srcsize;++j) {
              const value_type *__restrict__ ss = src + j*srcsize;
              /* */ value_type *__restrict__ tt = tgt + j*srcsize;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
              for (ordinal_type i=0;i<(j+1);++i)
                tt[i] += ss[i];
            }

            // unlock
            s.lock = 0;
            Kokkos::load_fence();
              
            return 0;
          }
        } 
        
        const size_type s2tsize = srcsize*sizeof(ordinal_type);
        TACHO_TEST_FOR_ABORT(bufsize < s2tsize, "bufsize is smaller than required s2t workspace");        

        dense_block_type A;
        ordinal_type *s2t = (ordinal_type*)buf;

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
            A.set_view(s.m, s.n);
            A.attach_buffer(1, s.m, s.buf);
            
            ordinal_type ijbeg = 0; for (;s2t[ijbeg] == -1; ++ijbeg) ;

            // lock
            while (Kokkos::atomic_compare_exchange(&s.lock, 0, 1)) KOKKOS_IMPL_PAUSE;
            Kokkos::store_fence();            

            for (ordinal_type jj=ijbeg;jj<srcsize;++jj) 
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
              for (ordinal_type ii=ijbeg;ii<srcsize;++ii) {
                const ordinal_type row = s2t[ii];
                if (row < s.m) A(row, s2t[jj]) += ABR(ii, jj);
                else break;
              }
            
            // unlock
            s.lock = 0;
            Kokkos::load_fence();
          }
        }
// #if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
// #if defined( TACHO_PROFILE_TIME_PER_THREAD )
//           g_time_per_thread[omp_get_thread_num()] += timer.seconds();
// #endif
// #endif
        return 0;
      }

      template<typename SchedulerType,
               typename MemberType,
               typename SupernodeInfoType>
      KOKKOS_INLINE_FUNCTION
      static int
      solve_lower(const SchedulerType &sched,
                  const MemberType &member,
                  const SupernodeInfoType &info,
                  const typename SupernodeInfoType::value_type_matrix &xB,
                  const ordinal_type sid) {
        typedef SupernodeInfoType supernode_info_type;

        typedef typename supernode_info_type::value_type value_type;
        typedef typename supernode_info_type::value_type_matrix value_type_matrix;

        typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

        const auto &s = info.supernodes(sid);

        // get panel pointer
        value_type *ptr = s.buf; 

        // panel is divided into diagonal and interface block
        const ordinal_type m = s.m, n = s.n - s.m, nrhs = info.x.dimension_1();

        // m and n are available, then factorize the supernode block
        if (m > 0) {
          const ordinal_type offm = s.row_begin;
          UnmanagedViewType<value_type_matrix> AL(ptr, m, m); ptr += m*m;
          auto xT = Kokkos::subview(info.x, range_type(offm, offm+m), Kokkos::ALL());

          if (nrhs >= ThresholdSolvePhaseUsingBlas3)
            Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::External>
              ::invoke(sched, member, Diag::NonUnit(), 1.0, AL, xT);
          else
            Trsv<Uplo::Upper,Trans::ConjTranspose,Algo::External>
              ::invoke(sched, member, Diag::NonUnit(), AL, xT);
            
          if (n > 0) {
            UnmanagedViewType<value_type_matrix> AR(ptr, m, n); // ptr += m*n;
            if (nrhs >= ThresholdSolvePhaseUsingBlas3)
              Gemm<Trans::ConjTranspose,Trans::NoTranspose,Algo::External>
                ::invoke(sched, member, -1.0, AR, xT, 0.0, xB);
            else
              Gemv<Trans::ConjTranspose,Algo::External>
                ::invoke(sched, member, -1.0, AR, xT, 0.0, xB);
          }
        }
        return 0;
      }

      template<typename SchedulerType,
               typename MemberType,
               typename SupernodeInfoType>
      KOKKOS_INLINE_FUNCTION
      static int
      update_solve_lower(const SchedulerType &sched,
                         const MemberType &member,
                         const SupernodeInfoType &info,
                         const typename SupernodeInfoType::value_type_matrix &xB,
                         const ordinal_type sid) {
        //typedef SupernodeInfoType supernode_info_type;
        //typedef typename supernode_info_type::value_type_matrix value_type_matrix;

        const auto &cur = info.supernodes(sid);
        const ordinal_type 
          sbeg = cur.sid_col_begin + 1, send = cur.sid_col_end - 1;

        const ordinal_type m = xB.dimension_0(), n = xB.dimension_1();
        TACHO_TEST_FOR_ABORT(m != (cur.n-cur.m), "# of rows in xB does not match to super blocksize in sid");
        
        for (ordinal_type i=sbeg,is=0;i<send;++i) {
          const ordinal_type 
            tbeg = info.sid_block_colidx(i).second,
            tend = info.sid_block_colidx(i+1).second;
          
          // lock
          const auto &s = info.supernodes(info.sid_block_colidx(i).first);
          while (Kokkos::atomic_compare_exchange(&s.lock, 0, 1)) KOKKOS_IMPL_PAUSE;
          Kokkos::store_fence();            
          
          // both src and tgt increase index
          for (ordinal_type it=tbeg;it<tend;++it,++is) {
            const ordinal_type row = info.gid_colidx(cur.gid_col_begin + it);
            for (ordinal_type j=0;j<n;++j) 
              info.x(row,j) += xB(is,j);
          }
          
          // unlock
          s.lock = 0;
          Kokkos::load_fence();          
        }
        
        return 0;
      }
      
      template<typename SchedulerType,
               typename MemberType,
               typename SupernodeInfoType>
      KOKKOS_INLINE_FUNCTION
      static int
      solve_upper(const SchedulerType &sched,
                  const MemberType &member,
                  const SupernodeInfoType &info,
                  const typename SupernodeInfoType::value_type_matrix &xB,
                  const ordinal_type sid) {
        typedef SupernodeInfoType supernode_info_type;

        typedef typename supernode_info_type::value_type value_type;
        typedef typename supernode_info_type::value_type_matrix value_type_matrix;

        typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

        // get current supernode
        const auto &s = info.supernodes(sid);

        // get supernode panel pointer
        value_type *ptr = s.buf;

        // panel is divided into diagonal and interface block
        const ordinal_type m = s.m, n = s.n - s.m, nrhs = info.x.dimension_1();

        // m and n are available, then factorize the supernode block
        if (m > 0) {
          const UnmanagedViewType<value_type_matrix> AL(ptr, m, m); ptr += m*m;

          const ordinal_type offm = s.row_begin;
          const auto xT = Kokkos::subview(info.x, range_type(offm, offm+m), Kokkos::ALL());

          if (n > 0) {
            const UnmanagedViewType<value_type_matrix> AR(ptr, m, n); // ptr += m*n;
            if (nrhs >= ThresholdSolvePhaseUsingBlas3)
              Gemm<Trans::NoTranspose,Trans::NoTranspose,Algo::External>
                ::invoke(sched, member, -1.0, AR, xB, 1.0, xT);
            else
              Gemv<Trans::NoTranspose,Algo::External>
                ::invoke(sched, member, -1.0, AR, xB, 1.0, xT);
          }
          if (nrhs >= ThresholdSolvePhaseUsingBlas3)
            Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,Algo::External>
              ::invoke(sched, member, Diag::NonUnit(), 1.0, AL, xT);
          else
            Trsv<Uplo::Upper,Trans::NoTranspose,Algo::External>
              ::invoke(sched, member, Diag::NonUnit(), AL, xT);
        }
        return 0;
      }

      template<typename SchedulerType,
               typename MemberType,
               typename SupernodeInfoType>
      KOKKOS_INLINE_FUNCTION
      static int
      update_solve_upper(const SchedulerType &sched,
                         const MemberType &member,
                         const SupernodeInfoType &info,
                         const typename SupernodeInfoType::value_type_matrix &xB,
                         const ordinal_type sid) {

        //typedef SupernodeInfoType supernode_info_type;
        //typedef typename supernode_info_type::value_type_matrix value_type_matrix;

        const auto &s = info.supernodes(sid);

        const ordinal_type m = xB.dimension_0(), n = xB.dimension_1();
        TACHO_TEST_FOR_ABORT(m != (s.n-s.m), "# of rows in xB does not match to super blocksize in sid");

        const ordinal_type goffset = s.gid_col_begin + s.m;
        for (ordinal_type j=0;j<n;++j)
          for (ordinal_type i=0;i<m;++i) {
            const ordinal_type row = info.gid_colidx(i+goffset);
            xB(i,j) = info.x(row,j);
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
                                 const ordinal_type bufsize) {
        typedef SupernodeInfoType supernode_info_type;

        typedef typename supernode_info_type::value_type value_type;
        typedef typename supernode_info_type::value_type_matrix value_type_matrix;

        const auto &s = info.supernodes(sid);

        if (final) {
          // serial recursion
          for (ordinal_type i=0;i<s.nchildren;++i)
            factorize_recursive_serial(sched, member, info, 
                                       s.children[i], final, buf, bufsize);
        }

        {
          const size_type n = s.n - s.m, bufsize_required = n*(n+1)*sizeof(value_type);

          TACHO_TEST_FOR_ABORT(bufsize < static_cast<ordinal_type>(bufsize_required), 
                               "bufsize is smaller than required");

          UnmanagedViewType<value_type_matrix> ABR((value_type*)buf, n, n);

          CholSupernodes<Algo::Workflow::Serial>
            ::factorize(sched, member, info, ABR, sid);

          CholSupernodes<Algo::Workflow::Serial>
            ::update(sched, member, info, ABR, sid,
                     bufsize - ABR.span()*sizeof(value_type),
                     (void*)((value_type*)buf + ABR.span()));
        }
        return 0;
      }


      template<typename SchedulerType,
               typename MemberType,
               typename SupernodeInfoType>
      KOKKOS_INLINE_FUNCTION
      static int
      solve_lower_recursive_serial(const SchedulerType &sched,
                                   const MemberType &member,
                                   const SupernodeInfoType &info,
                                   const ordinal_type sid,
                                   const bool final,
                                   typename SupernodeInfoType::value_type *buf,
                                   const ordinal_type bufsize) {
        typedef SupernodeInfoType supernode_info_type;
        
        typedef typename supernode_info_type::value_type value_type;
        typedef typename supernode_info_type::value_type_matrix value_type_matrix;
        
        const auto &s = info.supernodes(sid);
        
        if (final) {
          // serial recursion
          for (ordinal_type i=0;i<s.nchildren;++i)
            solve_lower_recursive_serial(sched, member, info, 
                                         s.children[i], final, buf, bufsize);
        }

        {
          const size_type n = s.n - s.m, nrhs = info.x.dimension_1(), 
            bufsize_required = n*nrhs*sizeof(value_type);

          TACHO_TEST_FOR_ABORT(bufsize < static_cast<ordinal_type>(bufsize_required), 
                               "bufsize is smaller than required");

          UnmanagedViewType<value_type_matrix> xB((value_type*)buf, n, nrhs);

          CholSupernodes<Algo::Workflow::Serial>
            ::solve_lower(sched, member, info, xB, sid);

          CholSupernodes<Algo::Workflow::Serial>
            ::update_solve_lower(sched, member, info, xB, sid);
        }
        return 0;
      }


      template<typename SchedulerType,
               typename MemberType,
               typename SupernodeInfoType>
      KOKKOS_INLINE_FUNCTION
      static int
      solve_upper_recursive_serial(const SchedulerType &sched,
                                   const MemberType &member,
                                   const SupernodeInfoType &info,
                                   const ordinal_type sid,
                                   const bool final,
                                   typename SupernodeInfoType::value_type *buf,
                                   const ordinal_type bufsize) {
        typedef SupernodeInfoType supernode_info_type;
        
        typedef typename supernode_info_type::value_type value_type;
        typedef typename supernode_info_type::value_type_matrix value_type_matrix;


        const auto &s = info.supernodes(sid);        
        {
          const size_type n = s.n - s.m, nrhs = info.x.dimension_1(), 
            bufsize_required = n*nrhs*sizeof(value_type);
          
          TACHO_TEST_FOR_ABORT(bufsize < static_cast<ordinal_type>(bufsize_required), 
                               "bufsize is smaller than required");

          UnmanagedViewType<value_type_matrix> xB((value_type*)buf, n, nrhs);

          CholSupernodes<Algo::Workflow::Serial>
            ::update_solve_upper(sched, member, info, xB, sid);

          CholSupernodes<Algo::Workflow::Serial>
            ::solve_upper(sched, member, info, xB, sid);
        }

        if (final) {
          // serial recursion
          for (ordinal_type i=0;i<s.nchildren;++i)
            solve_upper_recursive_serial(sched, member, info, 
                                         s.children[i], final, buf, bufsize);
        }
        return 0;
      }

    };
  }
}

#endif
