#ifndef __TACHO_SUPERNODE_INFO_HPP__
#define __TACHO_SUPERNODE_INFO_HPP__

#include "Tacho_Util.hpp"

/// \file Tacho_SupernodeInfo.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

    struct SuperNodeInfoInitReducer {
      typedef SuperNodeInfoInitReducer reducer;
      struct ValueType {
        ordinal_type max_nchildren;
        ordinal_type max_supernode_size;
        ordinal_type max_schur_size;
        size_type nnz;
      };
      typedef struct ValueType value_type;
        
      typedef Kokkos::View<value_type,Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;
        
      value_type *value;
        
      KOKKOS_INLINE_FUNCTION SuperNodeInfoInitReducer() = default;
      KOKKOS_INLINE_FUNCTION SuperNodeInfoInitReducer(const SuperNodeInfoInitReducer &b) = default;
      KOKKOS_INLINE_FUNCTION SuperNodeInfoInitReducer(value_type &val) : value(&val) {}
        
      KOKKOS_INLINE_FUNCTION void join(value_type &dst, value_type &src) const {
        dst.max_nchildren = ( src.max_nchildren > dst.max_nchildren ?
                              src.max_nchildren : dst.max_nchildren );
        dst.max_supernode_size = ( src.max_supernode_size > dst.max_supernode_size ?
                                   src.max_supernode_size : dst.max_supernode_size );
        dst.max_schur_size = ( src.max_schur_size > dst.max_schur_size ?
                               src.max_schur_size : dst.max_schur_size );
        dst.nnz += src.nnz;
      }
        
      KOKKOS_INLINE_FUNCTION void join(volatile value_type &dst, const volatile value_type &src) const {
        dst.max_nchildren = ( src.max_nchildren > dst.max_nchildren ?
                              src.max_nchildren : dst.max_nchildren );
        dst.max_supernode_size = ( src.max_supernode_size > dst.max_supernode_size ?
                                   src.max_supernode_size : dst.max_supernode_size );
        dst.max_schur_size = ( src.max_schur_size > dst.max_schur_size ?
                               src.max_schur_size : dst.max_schur_size );
        dst.nnz += src.nnz;
      }
        
      KOKKOS_INLINE_FUNCTION void init(value_type &val) const {
        val.max_nchildren = Kokkos::reduction_identity<ordinal_type>::max();
        val.max_supernode_size = Kokkos::reduction_identity<ordinal_type>::max();
        val.max_schur_size = Kokkos::reduction_identity<ordinal_type>::max();
        val.nnz = Kokkos::reduction_identity<size_type>::sum();
      }
        
      KOKKOS_INLINE_FUNCTION
      value_type& reference() {
        return *value;
      }
        
      KOKKOS_INLINE_FUNCTION
      result_view_type view() const {
        return result_view_type(value);
      }
    };

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
        ordinal_type nchildren, *children; //children[MaxDependenceSize]; // hierarchy

        ordinal_type max_decendant_schur_size;      // workspace
        ordinal_type max_decendant_supernode_size;  // workspace

        value_type *buf;

        KOKKOS_INLINE_FUNCTION
        Supernode() 
          : lock(0), row_begin(0), m(0), n(0), 
            gid_col_begin(0), gid_col_end(0), 
            sid_col_begin(0), sid_col_end(0), 
            nchildren(0), children(NULL),
            max_decendant_schur_size(0),
            max_decendant_supernode_size(0),
            buf(NULL) {
          //for (ordinal_type i=0;i<MaxDependenceSize;++i) children[i] = 0;
        }

        KOKKOS_INLINE_FUNCTION
        Supernode(const Supernode &b) 
          : lock(0), row_begin(b.row_begin), m(b.m), n(b.n), 
            gid_col_begin(b.gid_col_begin), gid_col_end(b.gid_col_end), 
            sid_col_begin(b.sid_col_begin), sid_col_end(b.sid_col_end), 
            nchildren(b.nchildren), children(b.children),
            max_decendant_schur_size(b.max_decendant_schur_size),
            max_decendant_supernode_size(b.max_decendant_supernode_size),
            buf(b.buf) {
          //for (ordinal_type i=0;i<b.nchildren;++i) children[i] = b.children[i];
        }

        KOKKOS_INLINE_FUNCTION
        ~Supernode()  {}        
      };
      typedef struct Supernode supernode_type;
      typedef Kokkos::View<supernode_type*,exec_space> supernode_type_array;

      ///
      /// info for symbolic
      ///
      //ConstUnmanagedViewType<supernode_type_array> supernodes;
      UnmanagedViewType<supernode_type_array> supernodes;

      /// dof mapping to sparse matrix
      UnmanagedViewType<ordinal_type_array> gid_colidx;

      /// supernode map and panel size configuration 
      /// first - sid, second - blk , blk_superpanel_colidx;
      /// the last sid is dummy but last blk is ending point of the block
      UnmanagedViewType<ordinal_pair_type_array> sid_block_colidx; 

      ///
      /// max parameter
      ///
      ordinal_type max_nchildren, max_supernode_size, max_schur_size;

      ///
      /// frontal matrix subassembly mode and serialization parameter
      ///
      short front_update_mode, serial_thres_size; // 0 - lock, 1 - atomic

      ///
      /// info for solve (rhs multivector)
      UnmanagedViewType<value_type_matrix> x;

      KOKKOS_INLINE_FUNCTION
      SupernodeInfo() = default;

      KOKKOS_INLINE_FUNCTION
      SupernodeInfo(const SupernodeInfo &b) = default;

      static 
      inline
      void
      initialize(/* */ SupernodeInfo &self,
                 /* */ supernode_type_array &supernodes_,
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
        const ordinal_type nsupernodes = snodes_.extent(0) - 1;

        /// allocate and assign supernodes
        supernodes_ = supernode_type_array("supernodes", nsupernodes); // managed view

        sid_block_colidx_ = ordinal_pair_type_array("sid_block_colidx", sid_colidx_.span());

        // by default, update mode is atomic: 0 - mutex lock, 1 - atomic
        self.front_update_mode = 1;
        self.serial_thres_size = 0;

        /// workspace parameter initialization
        self.max_nchildren = 0;
        self.max_supernode_size = 0;
        self.max_schur_size = 0; 
        
        Kokkos::RangePolicy<exec_space> supernodes_range_policy(0,nsupernodes);
        SuperNodeInfoInitReducer::value_type init_reduce_val;        
        Kokkos::parallel_reduce
          (supernodes_range_policy, 
           KOKKOS_LAMBDA(const ordinal_type &sid, SuperNodeInfoInitReducer::value_type &update) {
            auto &s = supernodes_(sid);
            
            s.row_begin = snodes_(sid);
            s.m = snodes_(sid+1) - snodes_(sid);
            s.n = blk_colidx_(sid_ptr_(sid+1)-1);
            
            s.gid_col_begin = gid_ptr_(sid); s.gid_col_end = gid_ptr_(sid+1);
            s.sid_col_begin = sid_ptr_(sid); s.sid_col_end = sid_ptr_(sid+1);
            
            for (ordinal_type i=s.sid_col_begin;i<s.sid_col_end;++i) {
              sid_block_colidx_(i).first  = sid_colidx_(i);
              sid_block_colidx_(i).second = blk_colidx_(i);
            }
            
            s.nchildren = stree_ptr_(sid+1) - stree_ptr_(sid);
            s.children = &stree_children_(stree_ptr_(sid));
            // const ordinal_type offset = stree_ptr_(sid);
            // for (ordinal_type i=0;i<s.nchildren;++i)
            //   s.children[i] = stree_children_(offset + i);

            update.max_nchildren = max(update.max_nchildren, s.nchildren);            
            update.max_supernode_size = max(update.max_supernode_size, s.m);
            update.max_schur_size = max(update.max_schur_size, s.n-s.m);
            update.nnz += s.m * s.n;
            
            s.max_decendant_supernode_size = s.m;
            s.max_decendant_schur_size = s.n-s.m;
            
          }, SuperNodeInfoInitReducer(init_reduce_val));

        {
          // finding max decendant supernode information is sequential
          // or level-based impl is required (too cumbersome..) 
          auto h_supernodes = Kokkos::create_mirror_view(supernodes_); 
          auto h_stree_parent = Kokkos::create_mirror_view(stree_parent_);           
          Kokkos::deep_copy(h_supernodes, supernodes_);
          Kokkos::deep_copy(h_stree_parent, stree_parent_);
          for (ordinal_type sid=0;sid<nsupernodes;++sid) {
            auto &s = h_supernodes(sid);
            const ordinal_type sidpar = h_stree_parent(sid);
            if (sidpar != -1) {
              auto &spar = h_supernodes(sidpar);
              spar.max_decendant_supernode_size = max(s.max_decendant_supernode_size,
                                                      spar.max_decendant_supernode_size);
              spar.max_decendant_schur_size = max(s.max_decendant_schur_size,
                                                  spar.max_decendant_schur_size);
            }
          }
          Kokkos::deep_copy(supernodes_, h_supernodes);
        }
        // need to iterate parallel for .. let's not do this way
        // Kokkos::parallel_for
        //   (supernodes_range_policy, KOKKOS_LAMBDA(const ordinal_type &sid) {
        //     auto &s = supernodes_(sid);
        //     const ordinal_type sidpar = stree_parent_(sid);
        //     if (sidpar != -1) {
        //       auto &spar = supernodes_(sidpar);
        //       spar.max_decendant_supernode_size = max(s.max_decendant_supernode_size,
        //                                               spar.max_decendant_supernode_size);
        //       spar.max_decendant_schur_size = max(s.max_decendant_schur_size,
        //                                           spar.max_decendant_schur_size);
        //     }
        //   });

        self.max_nchildren = init_reduce_val.max_nchildren;
        self.max_supernode_size = init_reduce_val.max_supernode_size;
        self.max_schur_size = init_reduce_val.max_schur_size;

        // supernodal factor array; data is held outside with a managed view
        // supernode does not include this view
        superpanel_buf_ = value_type_array("superpanel_buf", init_reduce_val.nnz); 
        Kokkos::parallel_scan
          (supernodes_range_policy, KOKKOS_LAMBDA(const ordinal_type &sid, size_type &update, const bool &final) {
            auto &s = supernodes_(sid);
            if (final) 
              s.buf = &superpanel_buf_(update);
            update += s.m * s.n;
          });

        self.supernodes  = supernodes_;  // unmanaged view, data is held outside
        self.gid_colidx = gid_colidx_;
        self.sid_block_colidx  = sid_block_colidx_;
      }

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
        initialize(*this,
                   supernodes_,
                   sid_block_colidx_,
                   superpanel_buf_,
                   snodes_,
                   gid_ptr_,
                   gid_colidx_,
                   sid_ptr_,
                   sid_colidx_,
                   blk_colidx_,
                   stree_parent_,
                   stree_ptr_,
                   stree_children_);
      }
      
      static 
      inline
      void
      copySparseToSuperpanels(SupernodeInfo &self, 
                              // input from sparse matrix
                              const size_type_array &ap,
                              const ordinal_type_array &aj,
                              const value_type_array &ax,
                              const ordinal_type_array &perm,
                              const ordinal_type_array &peri) {
#if 1
        const ordinal_type nsupernodes = self.supernodes.extent(0);
        Kokkos::TeamPolicy<exec_space,Kokkos::Schedule<Kokkos::Static> > 
          policy(nsupernodes, Kokkos::AUTO()); // team and vector sizes are AUTO selected.

        Kokkos::parallel_for
          (policy, KOKKOS_LAMBDA (const typename Kokkos::TeamPolicy<exec_space>::member_type &member) {
            const ordinal_type sid = member.league_rank();
            const auto s = self.supernodes(sid);
            dense_block_type tgt(s.buf, s.m, s.n);;            

            // row major access to sparse src
            Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, s.m), [&](const ordinal_type &i) {
                const ordinal_type 
                  ii = i + s.row_begin,  // row in U
                  row = perm(ii), kbeg = ap(row), kend = ap(row+1);   // row in A
                const ordinal_type kcnt = kend - kbeg;
                Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, kcnt), [&](const ordinal_type &kk) {
                    const ordinal_type k  = kk + kbeg;
                    const ordinal_type jj = peri(aj(k) /* col in A */); // col in U
                    if (ii <= jj) {
                      ordinal_type *first = self.gid_colidx.data() + s.gid_col_begin; 
                      ordinal_type *last  = self.gid_colidx.data() + s.gid_col_end;
                      ordinal_type *loc   = lower_bound(first, last, jj,
                                                        [](ordinal_type left, ordinal_type right) { 
                                                          return left < right; });                      
                      TACHO_TEST_FOR_ABORT(*loc != jj, " copy is wrong" );
                      tgt(i, loc-first) = ax(k);
                    }
                  });
              });
          });
#else
        const ordinal_type nsupernodes = self.supernodes.extent(0);
        const ordinal_type m = ap.extent(0) - 1;
        Kokkos::TeamPolicy<exec_space,Kokkos::Schedule<Kokkos::Static> > 
          policy(nsupernodes, Kokkos::AUTO()); // team and vector sizes are AUTO selected.

        typedef typename exec_space::scratch_memory_space shmem_space;        
        typedef Kokkos::View<ordinal_type*,shmem_space,Kokkos::MemoryUnmanaged> team_shared_memory_view_type;
        const ordinal_type lvl = 0, per_team_scratch = team_shared_memory_view_type::shmem_size(m);
        
        Kokkos::parallel_for
          (policy.set_scratch_size(lvl, Kokkos::PerTeam(per_team_scratch)),
           KOKKOS_LAMBDA ( const typename Kokkos::TeamPolicy<exec_space>::member_type &member) {
            team_shared_memory_view_type work(member.team_shmem(), m);
            const ordinal_type sid = member.league_rank();
            const auto s = self.supernodes(sid);
            dense_block_type tgt(s.buf, s.m, s.n);;            
            
            // local to global map
            Kokkos::parallel_for(Kokkos::TeamThreadRange(member, s.n), [&](const ordinal_type &j) {
                Kokkos::single(Kokkos::PerThread(member), [&]() {
                    work[self.gid_colidx(j+s.gid_col_begin) /* = col */] = j;
                  });
              });
            
            // row major access to sparse src
            Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, s.m), [&](const ordinal_type &i) {
                const ordinal_type 
                  ii = i + s.row_begin,  // row in U
                  row = perm(ii), kbeg = ap(row), kend = ap(row+1);   // row in A
                const ordinal_type kcnt = kend - kbeg;
                Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, kcnt), [&](const ordinal_type &kk) {
                    const ordinal_type k  = kk + kbeg;
                    const ordinal_type jj = peri(aj(k) /* col in A */); // col in U
                    if (ii <= jj) 
                      tgt(i, work[jj]) = ax(k);
                  });
              });
          });
#endif
      }

      inline
      void
      copySparseToSuperpanels(// input from sparse matrix
                              const size_type_array &ap,
                              const ordinal_type_array &aj,
                              const value_type_array &ax,
                              const ordinal_type_array &perm,
                              const ordinal_type_array &peri) {
        copySparseToSuperpanels(*this,
                                ap,
                                aj,
                                ax,
                                perm,
                                peri);      
      }      
      
      static
      inline
      void
      createCrsMatrix(SupernodeInfo &self,
                      crs_matrix_type &A,
                      const bool replace_value_with_one = false) {
        // count m, n, nnz
        const ordinal_type nsupernodes = self.supernodes.extent(0);

        auto d_last = Kokkos::subview(self.supernodes, nsupernodes - 1);
        auto h_last = Kokkos::create_mirror_view(d_last);
        Kokkos::deep_copy(h_last, d_last);
        auto last = h_last();
        //auto &last = supernodes(nsupernodes - 1);

        const ordinal_type mm = last.row_begin + last.m, nn = mm;

        Kokkos::RangePolicy<exec_space> supernodes_range_policy(0,nsupernodes);

        // parallel for/scan version
        size_type_array ap_tmp("ap", mm+1);
        Kokkos::parallel_for
          (supernodes_range_policy, KOKKOS_LAMBDA(const ordinal_type &sid) {
            // row major access to sparse src
            const auto &s = self.supernodes(sid);
            const ordinal_type soffset = s.row_begin;
            for (ordinal_type i=0;i<s.m;++i) {
              const ordinal_type row = i+soffset; // row in sparse matrix
              ap_tmp(row) = (s.n - i); // upper triangular only
            }
          });

        size_type_array ap("ap", mm+1);
        Kokkos::parallel_scan
          (Kokkos::RangePolicy<exec_space>(0,mm+1), 
           KOKKOS_LAMBDA(const ordinal_type &i, size_type &update, const bool &final) {
            if (final) 
              ap(i) = update;
            update += ap_tmp(i);
          });
        
        // fill the matrix
        auto d_nnz = Kokkos::subview(ap, mm);
        auto h_nnz = Kokkos::create_mirror_view(d_nnz);        
        Kokkos::deep_copy(h_nnz, d_nnz);

        const auto nnz = h_nnz();
        ordinal_type_array aj("aj", nnz);
        value_type_array ax("ax", nnz);

        Kokkos::parallel_for
          (supernodes_range_policy, KOKKOS_LAMBDA(const ordinal_type &sid) {
            const auto &s = self.supernodes(sid);
            
            dense_block_type src;
            src.set_view(s.m, s.n);
            src.attach_buffer(1, s.m, s.buf);
            
            // row major access to sparse src
            const ordinal_type
              soffset = s.row_begin,
              goffset = s.gid_col_begin;
            
            for (ordinal_type i=0;i<s.m;++i) {
              const size_type beg = ap(i+soffset);
              for (ordinal_type j=i,k=0;j<s.n;++j,++k) {
                const ordinal_type col = self.gid_colidx(j+goffset);
                aj(beg+k) = col;
                ax(beg+k) = replace_value_with_one ? 1.0 : src(i, j);
              }
            }
          });

        // set triple to crs matrix
        A.clear();
        A.setExternalMatrix(mm, nn, nnz, ap, aj, ax);
      }

      inline
      void
      createCrsMatrix(crs_matrix_type &A,
                      const bool replace_value_with_one = false) {
        createCrsMatrix(*this,
                        A, replace_value_with_one);
      }

    };
    
}

#endif
