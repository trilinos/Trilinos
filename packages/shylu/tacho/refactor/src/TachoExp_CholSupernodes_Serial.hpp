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

            TACHO_TEST_FOR_ABORT(ABR.dimension_0() != n ||
                                 ABR.dimension_1() != n,
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
        typedef SupernodeInfoType supernode_info_type;
        
        typedef typename supernode_info_type::ordinal_type_array ordinal_type_array;
        typedef typename supernode_info_type::dense_block_type dense_block_type;

        const auto &cur = info.supernodes(sid);

        const ordinal_type 
          sbeg = cur.sid_col_begin + 1, send = cur.sid_col_end - 1;

        const ordinal_type 
          srcbeg  = info.sid_block_colidx(sbeg).second, 
          srcend  = info.sid_block_colidx(send).second, 
          srcsize = srcend - srcbeg;

        const size_type s2tsize = srcsize*sizeof(ordinal_type);
        TACHO_TEST_FOR_ABORT(bufsize < s2tsize, "bufsize is smaller than required s2t workspace");        

        dense_block_type A;
        ordinal_type *s2t = (ordinal_type*)buf;

        // loop over target
        const ordinal_type *s_colidx = sbeg < send ? &info.gid_colidx(cur.gid_col_begin + srcbeg) : NULL;
        for (ordinal_type i=sbeg;i<send;++i) {
          const auto &s = info.supernodes(info.sid_block_colidx(i).first);

          A.set_view(s.m, s.n);
          A.attach_buffer(1, s.m, s.buf);

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

          {
            ordinal_type ijbeg = 0; for (;s2t[ijbeg] == -1; ++ijbeg) ;
            for (ordinal_type jj=ijbeg;jj<srcsize;++jj) 
              for (ordinal_type ii=ijbeg;ii<srcsize;++ii) {
                const ordinal_type row = s2t[ii];
                if (row < s.m) 
                  A(row, s2t[jj]) += ABR(ii, jj);
                else
                  break;
              }
          }

          // ordinal_type ijbeg = 0; for (;s2t[ijbeg] == -1; ++ijbeg) ;
          // for (ordinal_type ii=ijbeg;ii<srcsize;++ii) {
          //   const ordinal_type row = s2t[ii];
          //   if (row < s.m) 
          //     for (ordinal_type jj=ijbeg;jj<srcsize;++jj) ;
          //   A(row, s2t[jj]) += ABR(ii, jj);
          //   else
          //     break;
          // }

          // for (ordinal_type ii=0;ii<srcsize;++ii) {
          //   const ordinal_type row = s2t[ii];
          //   if (row < s.m && row != -1) {
          //     for (ordinal_type jj=0;jj<srcsize;++jj) {
          //       const ordinal_type col = s2t[jj];
          //       if (col != -1) 
          //         A(row, col) += ABR(ii, jj);
          //     }
          //   }
          // }
        }
        return 0;
      }

//       template<typename SchedulerType,
//                typename MemberType,
//                typename SupernodeInfoType>
//       KOKKOS_INLINE_FUNCTION
//       static int
//       solve_lower(const SchedulerType &sched,
//                   const MemberType &member,
//                   const SupernodeInfoType &info,
//                   const typename SupernodeInfoType::value_type_matrix &xB,
//                   const ordinal_type sid) {
//         typedef SupernodeInfoType supernode_info_type;

//         typedef typename supernode_info_type::value_type value_type;
//         typedef typename supernode_info_type::value_type_matrix value_type_matrix;

//         typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

//         // get supernode panel pointer
//         value_type *ptr = info.getSuperPanelPtr(sid);
//         const UnmanagedViewType<value_type_matrix> &x = info.x;

//         // characterize the panel size
//         ordinal_type pm, pn;
//         info.getSuperPanelSize(sid, pm, pn);

//         // panel is divided into diagonal and interface block
//         const ordinal_type m = pm, n = pn - pm;

//         // m and n are available, then factorize the supernode block
//         if (m > 0) {
//           const ordinal_type offm = info.supernodes(sid);
//           UnmanagedViewType<value_type_matrix> AL(ptr, m, m); ptr += m*m;
//           auto xT = Kokkos::subview(x, range_type(offm, offm+m), Kokkos::ALL());

//           Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::External>
//             ::invoke(sched, member, Diag::NonUnit(), 1.0, AL, xT);

//           if (n > 0) {
//             UnmanagedViewType<value_type_matrix> AR(ptr, m, n); // ptr += m*n;

//             TACHO_TEST_FOR_ABORT(n != xB.dimension_0(), "# of rows in xB does not match to super blocksize in sid");
//             TACHO_TEST_FOR_ABORT(xB.dimension_1() != x.dimension_1(),
//                                  "# of cols in xB does not match to rhs in supernodesinfo");

//             Gemm<Trans::ConjTranspose,Trans::NoTranspose,Algo::External>
//               ::invoke(sched, member, -1.0, AR, xT, 0.0, xB);
//           }
//         }
//         return 0;
//       }

//       template<typename SchedulerType,
//                typename MemberType,
//                typename SupernodeInfoType>
//       KOKKOS_INLINE_FUNCTION
//       static int
//       update_solve_lower(const SchedulerType &sched,
//                          const MemberType &member,
//                          const SupernodeInfoType &info,
//                          const typename SupernodeInfoType::value_type_matrix &xB,
//                          const ordinal_type sid) {
//         typedef SupernodeInfoType supernode_info_type;
//         typedef typename supernode_info_type::value_type_matrix value_type_matrix;

//         // grab super panel
//         ordinal_type pm, pn;
//         info.getSuperPanelSize(sid, pm, pn);

//         const ordinal_type m = xB.dimension_0(), n = xB.dimension_1();
//         TACHO_TEST_FOR_ABORT(m != (pn-pm), "# of rows in xB does not match to super blocksize in sid");

//         const UnmanagedViewType<value_type_matrix> &x = info.x;
//         TACHO_TEST_FOR_ABORT(n != x.dimension_1(), "# of cols in xB does not match to rhs in supernodesinfo");

//         const ordinal_type goffset = info.gid_super_panel_ptr(sid) + pm;
// #if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
// #pragma omp critical
//         {
//           for (ordinal_type j=0;j<n;++j)
//             for (ordinal_type i=0;i<m;++i) {
//               const ordinal_type row = info.gid_super_panel_colidx(i+goffset);
//               //Kokkos::atomic_fetch_add(&x(row, j), xB(i,j));
//               x(row,j) += xB(i,j);
//             }
//         }
// #else
//         {
//           for (ordinal_type j=0;j<n;++j)
//             for (ordinal_type i=0;i<m;++i) {
//               const ordinal_type row = info.gid_super_panel_colidx(i+goffset);
//               Kokkos::atomic_fetch_add(&x(row, j), xB(i,j));
//               //x(row,j) += xB(i,j);
//             }
//         }
// #endif

//         return 0;
//       }

//       template<typename SchedulerType,
//                typename MemberType,
//                typename SupernodeInfoType>
//       KOKKOS_INLINE_FUNCTION
//       static int
//       solve_upper(const SchedulerType &sched,
//                   const MemberType &member,
//                   const SupernodeInfoType &info,
//                   const typename SupernodeInfoType::value_type_matrix &xB,
//                   const ordinal_type sid) {
//         typedef SupernodeInfoType supernode_info_type;

//         typedef typename supernode_info_type::value_type value_type;
//         typedef typename supernode_info_type::value_type_matrix value_type_matrix;

//         typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

//         // get supernode panel pointer
//         value_type *ptr = info.getSuperPanelPtr(sid);
//         const UnmanagedViewType<value_type_matrix> &x = info.x;

//         // characterize the panel size
//         ordinal_type pm, pn;
//         info.getSuperPanelSize(sid, pm, pn);

//         // panel is divided into diagonal and interface block
//         const ordinal_type m = pm, n = pn - pm;

//         // m and n are available, then factorize the supernode block
//         if (m > 0) {
//           TACHO_TEST_FOR_ABORT(xB.dimension_1() != x.dimension_1(),
//                                "# of cols in xB does not match to rhs in supernodesinfo");

//           const UnmanagedViewType<value_type_matrix> AL(ptr, m, m); ptr += m*m;

//           const ordinal_type offm = info.supernodes(sid);
//           const auto xT = Kokkos::subview(x, range_type(offm, offm+m), Kokkos::ALL());

//           if (n > 0) {
//             const UnmanagedViewType<value_type_matrix> AR(ptr, m, n); // ptr += m*n;
//             Gemm<Trans::NoTranspose,Trans::NoTranspose,Algo::External>
//               ::invoke(sched, member, -1.0, AR, xB, 1.0, xT);
//           }
//           Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,Algo::External>
//             ::invoke(sched, member, Diag::NonUnit(), 1.0, AL, xT);
//         }
//         return 0;
//       }

//       template<typename SchedulerType,
//                typename MemberType,
//                typename SupernodeInfoType>
//       KOKKOS_INLINE_FUNCTION
//       static int
//       update_solve_upper(const SchedulerType &sched,
//                          const MemberType &member,
//                          const SupernodeInfoType &info,
//                          const typename SupernodeInfoType::value_type_matrix &xB,
//                          const ordinal_type sid) {

//         typedef SupernodeInfoType supernode_info_type;
//         typedef typename supernode_info_type::value_type_matrix value_type_matrix;

//         // grab super panel
//         ordinal_type pm, pn;
//         info.getSuperPanelSize(sid, pm, pn);

//         const ordinal_type m = xB.dimension_0(), n = xB.dimension_1();
//         TACHO_TEST_FOR_ABORT(m != (pn-pm), "# of rows in xB does not match to super blocksize in sid");

//         const UnmanagedViewType<value_type_matrix> &x = info.x;

//         const ordinal_type goffset = info.gid_super_panel_ptr(sid) + pm;
//         for (ordinal_type j=0;j<n;++j)
//           for (ordinal_type i=0;i<m;++i) {
//             const ordinal_type row = info.gid_super_panel_colidx(i+goffset);
//             xB(i,j) = x(row,j);
//           }

//         return 0;
//       }

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

          TACHO_TEST_FOR_ABORT(bufsize < bufsize_required, "bufsize is smaller than required");

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


      // template<typename SchedulerType,
      //          typename MemberType,
      //          typename SupernodeInfoType>
      // KOKKOS_INLINE_FUNCTION
      // static int
      // solve_lower_recursive_serial(const SchedulerType &sched,
      //                              const MemberType &member,
      //                              const SupernodeInfoType &info,
      //                              const ordinal_type sid,
      //                              const bool final,
      //                              typename SupernodeInfoType::value_type *buf,
      //                              const ordinal_type bufsize) {
      //   typedef SupernodeInfoType supernode_info_type;
        
      //   typedef typename supernode_info_type::value_type value_type;
      //   typedef typename supernode_info_type::value_type_matrix value_type_matrix;

      //   if (final) {
      //     // serial recursion
      //     const ordinal_type
      //       ibeg = info.stree_ptr(sid),
      //       iend = info.stree_ptr(sid+1);

      //     for (ordinal_type i=ibeg;i<iend;++i)
      //       solve_lower_recursive_serial(sched, member, info, info.stree_children(i), final, buf, bufsize);
      //   }

      //   {
      //     ordinal_type pm, pn; info.getSuperPanelSize(sid, pm, pn);                                             
      //     const size_type n = pn - pm, nrhs = info.x.dimension_1(), 
      //       bufsize_required = n*nrhs*sizeof(value_type);
          
      //     TACHO_TEST_FOR_ABORT(bufsize < bufsize_required, "bufsize is smaller than required");

      //     UnmanagedViewType<value_type_matrix> xB((value_type*)buf, n, nrhs);

      //     CholSupernodes<Algo::Workflow::Serial>
      //       ::solve_lower(sched, member, info, xB, sid);

      //     CholSupernodes<Algo::Workflow::Serial>
      //       ::update_solve_lower(sched, member, info, xB, sid);
      //   }
      //   return 0;
      // }


      // template<typename SchedulerType,
      //          typename MemberType,
      //          typename SupernodeInfoType>
      // KOKKOS_INLINE_FUNCTION
      // static int
      // solve_upper_recursive_serial(const SchedulerType &sched,
      //                              const MemberType &member,
      //                              const SupernodeInfoType &info,
      //                              const ordinal_type sid,
      //                              const bool final,
      //                              typename SupernodeInfoType::value_type *buf,
      //                              const ordinal_type bufsize) {
      //   typedef SupernodeInfoType supernode_info_type;
        
      //   typedef typename supernode_info_type::value_type value_type;
      //   typedef typename supernode_info_type::value_type_matrix value_type_matrix;

      //   {
      //     ordinal_type pm, pn; info.getSuperPanelSize(sid, pm, pn);                                             
      //     const size_type n = pn - pm, nrhs = info.x.dimension_1(), 
      //       bufsize_required = n*nrhs*sizeof(value_type);
          
      //     TACHO_TEST_FOR_ABORT(bufsize < bufsize_required, "bufsize is smaller than required");

      //     UnmanagedViewType<value_type_matrix> xB((value_type*)buf, n, nrhs);

      //     CholSupernodes<Algo::Workflow::Serial>
      //       ::update_solve_upper(sched, member, info, xB, sid);

      //     CholSupernodes<Algo::Workflow::Serial>
      //       ::solve_upper(sched, member, info, xB, sid);
      //   }

      //   if (final) {
      //     // serial recursion
      //     const ordinal_type
      //       ibeg = info.stree_ptr(sid),
      //       iend = info.stree_ptr(sid+1);

      //     for (ordinal_type i=ibeg;i<iend;++i)
      //       solve_upper_recursive_serial(sched, member, info, info.stree_children(i), final, buf, bufsize);
      //   }
      //   return 0;
      // }

    };
  }
}

#endif
