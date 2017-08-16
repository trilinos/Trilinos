#ifndef __TACHOEXP_SUPERNODE_INFO_HPP__
#define __TACHOEXP_SUPERNODE_INFO_HPP__

#include "TachoExp_Util.hpp"

/// \file TachoExp_SupernodeInfo.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  namespace Experimental {

    template<typename ValueType, typename ExecSpace>
    struct SupernodeInfo {
      typedef ValueType value_type;
      typedef ExecSpace exec_space;

      typedef CrsMatrixBase<value_type,exec_space> crs_matrix_type;

      typedef Kokkos::View<ordinal_type*,exec_space> ordinal_type_array;
      typedef Kokkos::View<size_type*,exec_space> size_type_array;
      typedef Kokkos::View<value_type*,exec_space> value_type_array;

      typedef Kokkos::pair<ordinal_type,ordinal_type> ordinal_pair_type;
      typedef Kokkos::View<ordinal_pair_type*,exec_space> ordinal_pair_type_array;
      typedef Kokkos::View<value_type**,Kokkos::LayoutLeft,exec_space> value_type_matrix;

      typedef DenseMatrixView<value_type,exec_space> dense_block_type;
      typedef DenseMatrixView<dense_block_type,exec_space> dense_matrix_of_blocks_type;
      
      typedef Kokkos::Future<int,exec_space> future_type;

      struct Supernode {
        mutable int32_t lock;

        ordinal_type row_begin;                     // beginning row
        ordinal_type m, n;                          // panel dimension

        // column connectivity (gid - dof, sid - supernode)
        ordinal_type gid_col_begin, gid_col_end, sid_col_begin, sid_col_end;  
        ordinal_type nchildren, children[MaxDependenceSize]; // hierarchy

        ordinal_type max_decendant_schur_size;      // workspace
        ordinal_type max_decendant_supernode_size;  // workspace

        value_type *buf;

        KOKKOS_INLINE_FUNCTION
        Supernode() 
          : lock(0), row_begin(0), m(0), n(0), 
            gid_col_begin(0), sid_col_begin(0), nchildren(0),
            max_decendant_schur_size(0),
            max_decendant_supernode_size(0),
            buf(NULL) {}

        KOKKOS_INLINE_FUNCTION
        Supernode(const Supernode &b) 
          : lock(0), row_begin(b.row_begin), m(b.m), n(b.n), 
            gid_col_begin(b.gid_col_begin), sid_col_begin(b.sid_col_begin), nchildren(b.nchildren),
            max_decendant_schur_size(b.max_decendant_schur_size),
            max_decendant_supernode_size(b.max_decendant_supernode_size),
            buf(b.buf) {
          for (ordinal_type i=0;i<b.nchildren;++i) children[i] = b.children[i];
        }
      };
      typedef struct Supernode supernode_type;
      typedef Kokkos::View<supernode_type*,exec_space> supernode_type_array;

      ///
      /// Phase 1: symbolic
      ///
      ConstUnmanagedViewType<supernode_type_array> supernodes;

      /// dof mapping to sparse matrix
      ConstUnmanagedViewType<ordinal_type_array> gid_colidx;

      /// supernode map and panel size configuration 
      /// first - sid, second - blk , blk_superpanel_colidx;
      /// the last sid is dummy but last blk is ending point of the block
      ConstUnmanagedViewType<ordinal_pair_type_array> sid_block_colidx; 

      ///
      /// Phase 2: max parameter
      ///
      ordinal_type max_supernode_size, max_schur_size, serial_thres_size;

      ///
      /// Phase 3: solve (rhs multivector)
      UnmanagedViewType<value_type_matrix> x;

      KOKKOS_INLINE_FUNCTION
      SupernodeInfo() = default;

      KOKKOS_INLINE_FUNCTION
      SupernodeInfo(const SupernodeInfo &b) = default;

      inline
      void
      initialize(/* */ supernode_type_array &supernodes_,
                 /* */ ordinal_pair_type_array &sid_block_colidx_,
                 /* */ value_type_array &superpanel_buf_,
                 // symbolic input
                 const ordinal_type_array &snodes_,
                 const size_type_array &gid_ptr_,
                 const ordinal_type_array &gid_colidx_,
                 const size_type_array &sid_ptr_,
                 const ordinal_type_array &sid_colidx_,
                 const ordinal_type_array &blk_colidx_,
                 // tree hierarchy
                 const ordinal_type_array &stree_parent_,
                 const size_type_array &stree_ptr_,
                 const ordinal_type_array &stree_children_) {
        const ordinal_type nsupernodes = snodes_.dimension_0() - 1;

        /// allocate and assign supernodes
        supernodes_ = supernode_type_array("supernodes", nsupernodes); // managed view
        supernodes  = supernodes_;  // unmanaged view, data is held outside

        gid_colidx        = gid_colidx_;

        sid_block_colidx_ = ordinal_pair_type_array("sid_block_colidx", sid_colidx_.span());
        sid_block_colidx  = sid_block_colidx_;

        /// workspace parameter initialization
        max_schur_size = 0; 
        max_supernode_size = 0;
        serial_thres_size = 0;

        size_type nnz = 0;
        for (ordinal_type sid=0;sid<nsupernodes;++sid) { 
          auto &s = supernodes_(sid);

          s.row_begin = snodes_(sid);
          s.m = snodes_(sid+1) - snodes_(sid);
          s.n = blk_colidx_(sid_ptr_(sid+1)-1);

          s.gid_col_begin = gid_ptr_(sid); s.gid_col_end = gid_ptr_(sid+1);
          s.sid_col_begin = sid_ptr_(sid); s.sid_col_end = sid_ptr_(sid+1);

          for (ordinal_type i=s.sid_col_begin;i<static_cast<ordinal_type>(sid_ptr_(sid+1));++i) {
            sid_block_colidx_(i).first  = sid_colidx_(i);
            sid_block_colidx_(i).second = blk_colidx_(i);
          }

          s.nchildren = stree_ptr_(sid+1) - stree_ptr_(sid);
          TACHO_TEST_FOR_EXCEPTION(s.nchildren > MaxDependenceSize, std::runtime_error,
                                   "# of children is greater than max dependence (children) size");

          const ordinal_type offset = stree_ptr_(sid);
          for (ordinal_type i=0;i<s.nchildren;++i) 
            s.children[i] = stree_children_(offset + i);

          max_supernode_size = max(max_supernode_size, s.m);
          max_schur_size = max(max_schur_size, s.n-s.m);

          s.max_decendant_supernode_size = max(s.m, s.max_decendant_supernode_size);
          s.max_decendant_schur_size = max(s.n-s.m, s.max_decendant_schur_size);

          const ordinal_type sidpar = stree_parent_(sid);
          if (sidpar != -1) {
            auto &spar = supernodes_(sidpar);
            spar.max_decendant_supernode_size = max(s.max_decendant_supernode_size,  
                                                    spar.max_decendant_supernode_size);
            spar.max_decendant_schur_size = max(s.max_decendant_schur_size,  
                                                spar.max_decendant_schur_size);
          }
          nnz += s.m * s.n;
        }
        
        // by default, serialization is made when supernode is smaller than max/4
        //serial_thres_size = max(max_supernode_size, max_schur_size)/4;

        // supernodal factor array; data is held outside with a managed view 
        // supernode does not include this view
        superpanel_buf_ = value_type_array("superpanel_buf", nnz); nnz = 0;
        for (ordinal_type sid=0;sid<nsupernodes;++sid) {
          auto &s = supernodes_(sid);
          s.buf = &superpanel_buf_(nnz);
          nnz += s.m * s.n;
        }
      }
      
      inline
      void
      copySparseToSuperpanels(// input from sparse matrix
                              const size_type_array &ap,
                              const ordinal_type_array &aj,
                              const value_type_array &ax,
                              const ordinal_type_array &perm,
                              const ordinal_type_array &peri,
                              // work array to store map
                              Kokkos::MemoryPool<exec_space> &pool) {

        const ordinal_type nsupernodes = supernodes.dimension_0(), m = ap.dimension_0() - 1;
        Kokkos::TeamPolicy<exec_space,
          Kokkos::Schedule<Kokkos::Static> > policy(nsupernodes, 1 /* teamsize */, 1024 /* worksetsize */);
        
        const int lvl = 0, per_team_scratch = m*sizeof(ordinal_type);
        
        Kokkos::parallel_for
          (policy.set_scratch_size(lvl, Kokkos::PerTeam(per_team_scratch)),
           KOKKOS_LAMBDA ( const typename Kokkos::TeamPolicy<exec_space>::member_type &member) {
            typedef typename exec_space::scratch_memory_space shmem_space;
            Kokkos::View<ordinal_type*,shmem_space,Kokkos::MemoryUnmanaged> work(member.team_shmem(), m);

            const auto &s = supernodes(member.league_rank());

            const dense_block_type tgt(s.buf, s.m, s.n);;            
            
            // local to global map
            for (ordinal_type j=0;j<s.n;++j) 
              work[gid_colidx(j+s.gid_col_begin) /* = col */] = j;
            
            // row major access to sparse src
            for (ordinal_type i=0;i<s.m;++i) {
              const ordinal_type 
                ii = i + s.row_begin,  // row in U
                row = perm(ii), kbeg = ap(row), kend = ap(row+1);   // row in A
              for (ordinal_type k=kbeg;k<kend;++k) {
                const ordinal_type jj = peri(aj(k) /* col in A */); // col in U
                if (ii <= jj) 
                  tgt(i, work[jj]) = ax(k);
              }
            }
          });

        // const Kokkos::RangePolicy<exec_space,Kokkos::Schedule<Kokkos::Static> > policy(0, nsupernodes);
        // Kokkos::parallel_for                                                                                
        //   (policy, KOKKOS_LAMBDA(const int sid) {                                                                     
        //     const auto &s = supernodes(sid);

        //     const dense_block_type tgt(s.buf, s.m, s.n);;            
            
        //     ordinal_type *work = (ordinal_type*)pool.allocate(m*sizeof(ordinal_type));
        //     TACHO_TEST_FOR_ABORT(work == NULL, "memory pool allocation fails");

        //     // local to global map
        //     for (ordinal_type j=0;j<s.n;++j) 
        //       work[gid_colidx(j+s.gid_col_begin) /* = col */] = j;
            
        //     // row major access to sparse src
        //     for (ordinal_type i=0;i<s.m;++i) {
        //       const ordinal_type 
        //         ii = i + s.row_begin,  // row in U
        //         row = perm(ii), kbeg = ap(row), kend = ap(row+1);   // row in A
        //       for (ordinal_type k=kbeg;k<kend;++k) {
        //         const ordinal_type jj = peri(aj(k) /* col in A */); // col in U
        //         if (ii <= jj) 
        //           tgt(i, work[jj]) = ax(k);
        //       }
        //     }
        //     pool.deallocate(work, m*sizeof(ordinal_type));
        //   });    
      }
      
      inline
      crs_matrix_type
      createCrsMatrix(const bool replace_value_with_one = false) {
        // count m, n, nnz
        const ordinal_type nsupernodes = supernodes.dimension_0();
        auto &last = supernodes(nsupernodes - 1);
        
        const ordinal_type mm = last.row_begin + last.m, nn = mm;

        size_type cnt = 0;
        size_type_array ap("ap", mm+1);
        for (ordinal_type sid=0;sid<nsupernodes;++sid) {
          const auto &s = supernodes(sid);

          // row major access to sparse src
          const ordinal_type soffset = s.row_begin;
          for (ordinal_type i=0;i<s.m;++i) {
            const ordinal_type row = i+soffset; // row in sparse matrix
            ap(row) = cnt;
            cnt += (s.n - i); // upper triangular only
          }
        }
        ap(mm) = cnt;
        
        // fill the matrix
        const size_type nnz = cnt; cnt = 0;
        ordinal_type_array aj("aj", nnz);
        value_type_array ax("ax", nnz);

        dense_block_type src;
        for (ordinal_type sid=0;sid<nsupernodes;++sid) {
          const auto &s = supernodes(sid);

          src.set_view(s.m, s.n);
          src.attach_buffer(1, s.m, s.buf);

          // row major access to sparse src
          const ordinal_type
            soffset = s.row_begin,
            goffset = s.gid_col_begin;

          for (ordinal_type i=0;i<s.m;++i) {
            const size_type beg = ap(i+soffset);
            for (ordinal_type j=i,k=0;j<s.n;++j,++k) {
              const ordinal_type col = gid_colidx(j+goffset);
              aj(beg+k) = col;
              ax(beg+k) = replace_value_with_one ? 1.0 : src(i, j);
            }
          }
        }
        
        // set triple to crs matrix
        crs_matrix_type r_val;
        r_val.setExternalMatrix(mm, nn, nnz, ap, aj, ax);

        return r_val;
      }
    };
    
  }
}

#endif
