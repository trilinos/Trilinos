#ifndef __TACHOEXP_TASKFUNCTOR_CHOLESKY_SUPERNODES_HPP__
#define __TACHOEXP_TASKFUNCTOR_CHOLESKY_SUPERNODES_HPP__

/// \file TachoExp_TaskFunctor_FactorizeSuperNodes.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

namespace Tacho {

  namespace Experimental {

    template<typename MatValueType, typename ExecSpace>
    struct TaskFunctor_CholeskySuperNodes {
      typedef ExecSpace exec_space;

      // task scheduler/future
      typedef Kokkos::TaskScheduler<exec_space> sched_type;
      typedef typename sched_type::member_type member_type;
      typedef Kokkos::Future<int,exec_space> future_type;

      // memory pool
      typedef Kokkos::Experimental::MemoryPool<exec_space> memory_pool_type;

      // value types
      typedef int value_type; // functor return type
      typedef MatValueType mat_value_type; // matrix value type

      // supernodes
      typedef SuperNodeInfo<mat_value_type,exec_space> supernode_info_type;
      typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

    private:
      sched_type _sched;
      memory_pool_type _pool;

      supernode_info_type _info;
      ordinal_type _sid, _sidpar, _state;

      // KOKKOS_INLINE_FUNCTION
      // void
      // update(const policy_type &policy,
      //        const member_type &member,
      //        const memory_pool_type &pool,
      //        const ordinal_type sid,
      //        supernode_info_type::value_type_matrix &ABR) {
      //   const size_type
      //     sbeg = _info._sid_super_panel_ptr(sid)+1,
      //     send = _info._sid_super_panel_ptr(sid+1)-1;

      //   const ordinal_type
      //     src_col_beg = _info._blk_super_panel_colidx(sbeg),
      //     src_col_end = _info._blk_super_panel_colidx(send),
      //     src_col_size = src_col_end - src_col_beg;

      //   const size_type alloc_size = src_col_size*sizeof(ordinal_type);
      //   uintptr_t *ptr_alloced = pool.allocate(alloc_size);

      //   supernode_info_type::ordinal_type_array map(ptr_alloced, src_col_size);
      //   const ordinal_type smapoff = _info._gid_super_panel_ptr(sid);
      //   auto src_map = Kokkos::subview(_info._gid_super_panel_colidx,
      //                                  range_type(smapoff + src_col_beg,smapoff + src_col_end));

      //   // walk through source rows
      //   UnmanagedViewType<supernode_info_type::value_type_matrix> A;
      //   const ordinal_type src_row_offset = _info._blk_super_panel_colidx(sbeg);
      //   for (size_type i=sbeg;i<send;++i) {
      //     /// ** soruce rows
      //     const ordinal_type
      //       src_row_beg = _info._blk_super_panel_colidx(i),
      //       src_row_end = _info._blk_super_panel_colidx(i+1);

      //     /// ** target rows
      //     const ordinal_type row = _info._sid_super_panel_colidx(i);

      //     ordinal_type m, n;
      //     _info.getSuperPanelSize(row, m, n);
      //     _info.getSuperPanel(row, m, n, A);

      //     /// ** map
      //     const size_type
      //       rbeg = _info._sid_super_panel_ptr(row),
      //       rend = _info._sid_super_panel_ptr(row+1)-1;

      //     const ordinal_type
      //       tgt_col_beg = _info._blk_super_panel_colidx(rbeg),
      //       tgt_col_end = _info._blk_super_panel_colidx(rend),
      //       tgt_col_size = tgt_col_end - tgt_col_beg;

      //     const ordinal_type tmapoff = _info._gid_super_panel_ptr(row);
      //     auto tgt_map = Kokkos::subview(_info._gid_super_panel_colidx,
      //                                    range_type(tmapoff + tgt_col_beg, tmapoff + tgt_col_end));

      //     for (ordinal_type k=0,l=0;k<src_col_size;++k) {
      //       map(k) = -1;
      //       for (;l<tgt_col_size && tgt_map(l) <= src_map(k);++l)
      //         if (src_map(k) == tgt_map(l)) {
      //           map(k) = l;
      //           break;
      //         }
      //     }

      //     // release future _info._supernodes_future(row) done; task spawn
      //     ordinal_type mbeg = 0; for (;map(mbeg) == -1; ++mbeg) ;
      //     for (ordinal_type jj=mbeg;jj<src_col_size;++jj) {
      //       const ordinal_type mj = map(jj);
      //       for (ordinal_type ii=src_row_beg;ii<src_row_end;++ii) {
      //         const ordinal_type mi = map(ii-src_row_beg+mbeg);
      //         A(mi, mj) += ABR(ii-src_row_offset,jj);
      //       }
      //     }
      //   }
      //   pool.deallocate((void*)ptr_alloced, alloc_size);
      // }

      // KOKKOS_INLINE_FUNCTION
      // void
      // factorize(const policy_type &policy,
      //           const member_type &member,
      //           const memory_pool_type &pool,
      //           const ordinal_type sid,
      //           const ordinal_type sidpar) {
      //   value_type *ptr = _info.getSuperPanelPtr(sid);
      //   ordinal_type mm, nn;
      //   _info.getSuperPanelSize(sid, mm, nn);

      //   const ordinal_type m = mm, n = nn - mm;
      //   UnmanagedViewType<value_type_matrix> ATL(ptr, m, m), ATR(ptr+m*m, m, n);

      //   Chol<Uplo::Upper,Algo::External>
      //     ::invoke(policy, member, ATL);

      //   if (sidpar != -1) {
      //     Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::External>
      //       ::invoke(policy, member, Diag::NonUnit(), 1.0, ATL, ATR);

      //     const ordinal_type update_supernode_beg = _info._sid_super_panel_ptr(sid);
      //     const ordinal_type update_supernode_end = _info._sid_super_panel_ptr(sid+1);

      //     // diag, parent, empty one (range is constructed for block which needs end point)
      //     const bool is_direct_update = ((update_supernode_end - update_supernode_beg) == 3 &&
      //                                    _info._sid_super_panel_colidx(update_supernode_beg) == sidpar);
      //     if (is_direct_update) {
      //       UnmanagedViewType<value_type_matrix> ABR;
      //       _info.getSuperPanel(update_supernode_beg, n, n, ABR);
      //       Herk<Uplo::Upper,Trans::ConjTranspose,Algo::External>
      //         ::invoke(policy, member, -1.0, ATR, 1.0, ABR);
      //     } else {
      //       const size_type alloc_size = n*n*sizeof(typename supdernode_info_type::value_type);
      //       uintptr_t *ptr_alloced = pool.allocate(alloc_size);

      //       UnmanagedViewType<supernode_info_type::value_type_matrix> ABR(ptr_alloced, n, n);
      //       Herk<Uplo::Upper,Trans::ConjTranspose,Algo::External>
      //         ::invoke(policy, member, -1.0, ATR, 0.0, ABR);

      //       // copy back to its parent
      //       update(policy, member, pool, sid, ABR);

      //       pool.deallocate((void*)ptr_alloced, alloc_size);
      //     }
      //   }
      // }

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_CholeskySuperNodes() = delete;

      KOKKOS_INLINE_FUNCTION
      TaskFunctor_CholeskySuperNodes(const sched_type &sched,
                                     const memory_pool_type &pool,
                                     const supernode_info_type &info,
                                     const ordinal_type sid,
                                     const ordinal_type sidpar)                                     
        : _sched(sched),
          _pool(pool),
          _info(info),
          _sid(sid),
          _sidpar(sidpar),
          _state(0) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        if (get_team_rank(member) == 0) {

          bool is_execute = _state;
          if (_state == 0) {
            const ordinal_type 
              ibeg = _info.stree_ptr(_sid), 
              iend = _info.stree_ptr(_sid+1),
              isize = iend - ibeg;
            printf("sid = %d, sidpar = %d, nchild = %d, task is spawning its child\n", _sid, _sidpar, isize);
            if (isize > 0) {
              const size_type alloc_size = isize*sizeof(future_type);
              future_type *dep = (future_type*)_pool.allocate(alloc_size);
              
              for (ordinal_type i=0;i<isize;++i) {
                const ordinal_type child = _info.stree_children(i+ibeg);
                future_type f = Kokkos::task_spawn(Kokkos::TaskSingle(_sched, i ? Kokkos::TaskPriority::Low : Kokkos::TaskPriority::High),
                                                   TaskFunctor_CholeskySuperNodes(_sched, _pool, _info, 
                                                                                  child, _sid));
                dep[i] = f;
              }
              future_type dep_all;// = Kokkos::when_all(dep, isize);
              ++_state;
              printf("this task is respawned\n");
              Kokkos::respawn(this, dep_all, Kokkos::TaskPriority::High);
              _pool.deallocate((void*)dep, alloc_size);
            } else {
              is_execute = true;
            }
          }
          TACHO_TEST_FOR_ABORT(_state > 1, "deadlock : state is not finalized");
          if (is_execute) {
            printf("I would like to factorize this task here %d %d\n", _sid, _sidpar);
          }
        }
      }
    };
  }
}

#endif
