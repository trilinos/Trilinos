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
                const ordinal_type sid,
                const size_type bufsize,
                /* */ void *buf) {
        typedef SupernodeInfoType supernode_info_type;

        typedef typename supernode_info_type::value_type value_type;
        typedef typename supernode_info_type::value_type_matrix value_type_matrix;

        // get supernode panel pointer
        value_type *ptr = info.getSuperPanelPtr(sid);

        // characterize the panel size
        ordinal_type pm, pn;
        info.getSuperPanelSize(sid, pm, pn);

        // panel is divided into diagonal and interface block (i.e., ATL and ATR)
        const ordinal_type m = pm, n = pn - pm;

        // m and n are available, then factorize the supernode block
        if (m > 0) {
          UnmanagedViewType<value_type_matrix> ATL(ptr, m, m); ptr += m*m;
          Chol<Uplo::Upper,Algo::External>
            ::invoke(sched, member, ATL);

          if (n > 0) {
            UnmanagedViewType<value_type_matrix> ATR(ptr, m, n); // ptr += m*n;
            Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::External>
              ::invoke(sched, member, Diag::NonUnit(), 1.0, ATL, ATR);

            const size_type abrsize = n*n*sizeof(value_type);
            if (abrsize) {
              TACHO_TEST_FOR_ABORT(bufsize < abrsize, "bufsize is smaller than required schur workspace");
              UnmanagedViewType<value_type_matrix> ABR((value_type*)buf, n, n);

              Herk<Uplo::Upper,Trans::ConjTranspose,Algo::External>
                ::invoke(sched, member, -1.0, ATR, 0.0, ABR);
            }
          }
        }
        return 0;
      }

      template<typename SchedulerType,
               typename MemberType,
               typename SupernodeInfoType,
               typename MatrixViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      update(const SchedulerType &sched,
             const MemberType &member,
             const SupernodeInfoType &info,
             const MatrixViewType &ABR,
             const ordinal_type sid,
             const size_type bufsize,
             /* */ void *buf) {
        typedef SupernodeInfoType supernode_info_type;

        typedef typename supernode_info_type::value_type_matrix value_type_matrix;
        typedef typename supernode_info_type::ordinal_type_array ordinal_type_array;

        typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

        const size_type
          sbeg = info.sid_super_panel_ptr(sid)+1,
          send = info.sid_super_panel_ptr(sid+1)-1;

        const ordinal_type
          src_col_beg = info.blk_super_panel_colidx(sbeg),
          src_col_end = info.blk_super_panel_colidx(send),
          src_col_size = src_col_end - src_col_beg;

        const size_type mapsize = src_col_size*sizeof(ordinal_type);
        TACHO_TEST_FOR_ABORT(bufsize < mapsize, "bufsize is smaller than required map workspace");

        UnmanagedViewType<ordinal_type_array> map((ordinal_type*)buf, src_col_size);

        const ordinal_type smapoff = info.gid_super_panel_ptr(sid);
        auto src_map = Kokkos::subview(info.gid_super_panel_colidx,
                                       range_type(smapoff + src_col_beg,smapoff + src_col_end));

        // walk through source rows
        UnmanagedViewType<value_type_matrix> A;
        const ordinal_type src_row_offset = info.blk_super_panel_colidx(sbeg);
        for (size_type i=sbeg;i<send;++i) {
          /// ** soruce rows
          const ordinal_type
            src_row_beg = info.blk_super_panel_colidx(i),
            src_row_end = info.blk_super_panel_colidx(i+1);

          /// ** target rows
          const ordinal_type row = info.sid_super_panel_colidx(i);

          ordinal_type m, n;
          info.getSuperPanelSize(row, m, n);
          info.getSuperPanel(row, m, n, A);

          /// ** map
          const size_type
            rbeg = info.sid_super_panel_ptr(row),
            rend = info.sid_super_panel_ptr(row+1)-1;

          const ordinal_type
            tgt_col_beg = info.blk_super_panel_colidx(rbeg),
            tgt_col_end = info.blk_super_panel_colidx(rend),
            tgt_col_size = tgt_col_end - tgt_col_beg;

          const ordinal_type tmapoff = info.gid_super_panel_ptr(row);
          auto tgt_map = Kokkos::subview(info.gid_super_panel_colidx,
                                         range_type(tmapoff + tgt_col_beg, tmapoff + tgt_col_end));

          for (ordinal_type k=0,l=0;k<src_col_size;++k) {
            map(k) = -1;
            for (;l<tgt_col_size && tgt_map(l) <= src_map(k);++l)
              if (src_map(k) == tgt_map(l)) {
                map(k) = l;
                break;
              }
          }

          // release future info._supernodes_future(row) done; task spawn
          ordinal_type mbeg = 0; for (;map(mbeg) == -1; ++mbeg) ;
          for (ordinal_type jj=mbeg;jj<src_col_size;++jj) {
            const ordinal_type mj = map(jj);
            for (ordinal_type ii=src_row_beg;ii<src_row_end;++ii) {
              const ordinal_type mi = map(ii-src_row_beg+mbeg);

              // for parallel assembly, this should be used
              // TODO:: use row-wise atomic assembly later if this hurts performance on knl
              //        for gpu, this is the right way to do
              Kokkos::atomic_fetch_add(&A(mi, mj), ABR(ii-src_row_offset,jj));
              //A(mi, mj) += ABR(ii-src_row_offset,jj);
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
      solve_lower(const SchedulerType &sched,
                  const MemberType &member,
                  const SupernodeInfoType &info,
                  const typename SupernodeInfoType::value_type_matrix &xB,
                  const ordinal_type sid) {
        typedef SupernodeInfoType supernode_info_type;

        typedef typename supernode_info_type::value_type value_type;
        typedef typename supernode_info_type::value_type_matrix value_type_matrix;

        typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

        // get supernode panel pointer
        value_type *ptr = info.getSuperPanelPtr(sid);
        const UnmanagedViewType<value_type_matrix> &x = info.x;

        // characterize the panel size
        ordinal_type pm, pn;
        info.getSuperPanelSize(sid, pm, pn);

        // panel is divided into diagonal and interface block
        const ordinal_type m = pm, n = pn - pm;

        // m and n are available, then factorize the supernode block
        if (m > 0) {
          const ordinal_type offm = info.supernodes(sid);
          UnmanagedViewType<value_type_matrix> AL(ptr, m, m); ptr += m*m;
          auto xT = Kokkos::subview(x, range_type(offm, offm+m), Kokkos::ALL());

          Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::External>
            ::invoke(sched, member, Diag::NonUnit(), 1.0, AL, xT);

          if (n > 0) {
            UnmanagedViewType<value_type_matrix> AR(ptr, m, n); // ptr += m*n;

            TACHO_TEST_FOR_ABORT(n != xB.dimension_0(), "# of rows in xB does not match to super blocksize in sid");
            TACHO_TEST_FOR_ABORT(xB.dimension_1() != x.dimension_1(), 
                                 "# of cols in xB does not match to rhs in supernodesinfo");

            Gemm<Trans::ConjTranspose,Trans::NoTranspose,Algo::External>
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
        typedef SupernodeInfoType supernode_info_type;
        typedef typename supernode_info_type::value_type_matrix value_type_matrix;

        // grab super panel
        ordinal_type pm, pn;
        info.getSuperPanelSize(sid, pm, pn);

        const ordinal_type m = xB.dimension_0(), n = xB.dimension_1();
        TACHO_TEST_FOR_ABORT(m != (pn-pm), "# of rows in xB does not match to super blocksize in sid");

        const UnmanagedViewType<value_type_matrix> &x = info.x;
        TACHO_TEST_FOR_ABORT(n != x.dimension_1(), "# of cols in xB does not match to rhs in supernodesinfo");

        const ordinal_type goffset = info.gid_super_panel_ptr(sid) + pm;
        for (ordinal_type j=0;j<n;++j) 
          for (ordinal_type i=0;i<m;++i) {
            const ordinal_type row = info.gid_super_panel_colidx(i+goffset);
            Kokkos::atomic_fetch_add(&x(row, j), xB(i,j));
            //x(row,j) += xB(i,j);
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

        // get supernode panel pointer
        value_type *ptr = info.getSuperPanelPtr(sid);
        const UnmanagedViewType<value_type_matrix> &x = info.x;

        // characterize the panel size
        ordinal_type pm, pn;
        info.getSuperPanelSize(sid, pm, pn);

        // panel is divided into diagonal and interface block
        const ordinal_type m = pm, n = pn - pm;

        // m and n are available, then factorize the supernode block
        if (m > 0) {
          TACHO_TEST_FOR_ABORT(xB.dimension_1() != x.dimension_1(), 
                               "# of cols in xB does not match to rhs in supernodesinfo");

          const UnmanagedViewType<value_type_matrix> AL(ptr, m, m); ptr += m*m;

          const ordinal_type offm = info.supernodes(sid);
          const auto xT = Kokkos::subview(x, range_type(offm, offm+m), Kokkos::ALL());

          if (n > 0) {
            const UnmanagedViewType<value_type_matrix> AR(ptr, m, n); // ptr += m*n;
            Gemm<Trans::NoTranspose,Trans::NoTranspose,Algo::External>
              ::invoke(sched, member, -1.0, AR, xB, 1.0, xT);
          }
          Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,Algo::External>
            ::invoke(sched, member, Diag::NonUnit(), 1.0, AL, xT);
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

        typedef SupernodeInfoType supernode_info_type;
        typedef typename supernode_info_type::value_type_matrix value_type_matrix;
        
        // grab super panel
        ordinal_type pm, pn;
        info.getSuperPanelSize(sid, pm, pn);

        const ordinal_type m = xB.dimension_0(), n = xB.dimension_1();
        TACHO_TEST_FOR_ABORT(m != (pn-pm), "# of rows in xB does not match to super blocksize in sid");

        const UnmanagedViewType<value_type_matrix> &x = info.x;
        
        const ordinal_type goffset = info.gid_super_panel_ptr(sid) + pm;
        for (ordinal_type j=0;j<n;++j) 
          for (ordinal_type i=0;i<m;++i) {
            const ordinal_type row = info.gid_super_panel_colidx(i+goffset);
            xB(i,j) = x(row,j);
          }
        
        return 0;
      }

    };
  }
}

#endif
