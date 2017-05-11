#ifndef __TACHOEXP_CHOL_SUPERNODES_SERIAL_HPP__
#define __TACHOEXP_CHOL_SUPERNODES_SERIAL_HPP__

/// \file TachoExp_CholSuperNodes.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

namespace Tacho {

  namespace Experimental {

    template<>
    struct CholSuperNodes<Algo::Workflow::Serial> {
      template<typename SchedulerType,
               typename MemoryPoolType,
               typename MemberType,
               typename SuperNodeInfoType,
               typename MatrixViewType>
      KOKKOS_INLINE_FUNCTION
      static int 
      update(const SchedulerType &sched,
             const MemoryPoolType &pool, 
             const MemberType &member,
             const SuperNodeInfoType &info,
             const MatrixViewType &ABR,
             const ordinal_type sid,
             const size_type bufsize,
             /* */ void *buf) {
        typedef SuperNodeInfoType supernode_info_type;
        
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

        UnmanagedViewType<ordinal_type_array> map;
        const size_type mapsize = src_col_size*sizeof(ordinal_type);
        const ordinal_type diff = bufsize - mapsize;
        if (diff >= 0) {
          map = ordinal_type_array((ordinal_type*)buf, src_col_size);
        } else {
          TACHO_TEST_FOR_ABORT(true, "bufsize is smaller than requested");
          // ordinal_type *mapbuf = (ordinal_type*)pool.allocate(mapsize);
          // TACHO_TEST_FOR_ABORT(mapbuf == NULL, "pool allocation fails");
          // map = ordinal_type_array(mapbuf, src_col_size);
        }
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
              A(mi, mj) += ABR(ii-src_row_offset,jj);
            }
          }
        }

        if (diff >= 0) {
          // do nothing
        } else {
          //pool.deallocate((void*)map.data(), mapsize);
        }

        return 0;
      }

      template<typename SchedulerType,
               typename MemoryPoolType,
               typename MemberType,
               typename SuperNodeInfoType>
      KOKKOS_INLINE_FUNCTION
      static int 
      factorize(const SchedulerType &sched,
                const MemoryPoolType &pool, 
                const MemberType &member,
                const SuperNodeInfoType &info,
                const ordinal_type sid,
                const ordinal_type sidpar,
                const size_type bufsize,
                /* */ void *buf) {
        typedef SuperNodeInfoType supernode_info_type;

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
            ::invoke(sched, /* pool, */ member, ATL);
          
          if (n > 0) {
            UnmanagedViewType<value_type_matrix> ATR(ptr, m, n); // ptr += m*n;
            Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::External>
              ::invoke(sched, /* pool, */ member, Diag::NonUnit(), 1.0, ATL, ATR);
            
            const ordinal_type update_supernode_beg = info.sid_super_panel_ptr(sid);
            const ordinal_type update_supernode_end = info.sid_super_panel_ptr(sid+1);
            
            // diag, parent, (+1 because there is empty one; range is constructed for block which needs end point)
            const bool is_direct_update = ((update_supernode_end - update_supernode_beg) == 3 &&
                                           info.sid_super_panel_colidx(update_supernode_beg) == sidpar);
            UnmanagedViewType<value_type_matrix> ABR;
            if (is_direct_update) {
              ABR = value_type_matrix(info.getSuperPanelPtr(update_supernode_beg), n, n);
              Herk<Uplo::Upper,Trans::ConjTranspose,Algo::External>
                ::invoke(sched, /* pool, */ member, -1.0, ATR, 1.0, ABR);
            } else {
              const size_type abrsize = n*n*sizeof(value_type);
              const ordinal_type diff = bufsize - abrsize;
              if (diff > 0) {
                ABR = value_type_matrix((value_type*)buf, n, n);
              } else {
                TACHO_TEST_FOR_ABORT(true, "bufsize is smaller than requested");
                //value_type *abrbuf = (value_type*)pool.allocate(abrsize);
                //TACHO_TEST_FOR_ABORT(abrbuf == NULL, "pool allocation fails");
                //ABR = value_type_matrix(abrbuf, n, n);                
              }
              
              Herk<Uplo::Upper,Trans::ConjTranspose,Algo::External>
                ::invoke(sched, /* pool, */ member, -1.0, ATR, 0.0, ABR);

              update(sched, pool, member, info, ABR, sid, 
                     max(diff, 0), (diff > 0 ? (void*)((value_type*)buf + n*n) : NULL));
              
              if (diff > 0) {
                // do nothing
              } else {
                //pool.deallocate((void*)ABR.data(), abrsize);
              }
            }
          }
        }
        return 0;
      }
    };
  }
}

#endif
