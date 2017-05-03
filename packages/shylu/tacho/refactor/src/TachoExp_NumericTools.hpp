#ifndef __TACHOEXP_NUMERIC_TOOLS_HPP__
#define __TACHOEXP_NUMERIC_TOOLS_HPP__

/// \file TachoExp_NumericTools.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

#include "TachoExp_Chol.hpp"
#include "TachoExp_Chol_Upper_External.hpp"

#include "TachoExp_Trsm.hpp"
#include "TachoExp_Trsm_Left_Upper_ConjTrans_External.hpp"

#include "TachoExp_Herk.hpp"
#include "TachoExp_Herk_Upper_ConjTrans_External.hpp"

namespace Tacho {

  namespace Experimental {

    template<typename ValueType, typename ExecSpace>
    class NumericTools {
    public:
      typedef ValueType value_type;

      typedef Kokkos::DefaultHostExecutionSpace host_exec_space;
      typedef Kokkos::View<ordinal_type*,host_exec_space> ordinal_type_array_host;
      typedef Kokkos::View<size_type*,   host_exec_space> size_type_array_host;
      typedef Kokkos::View<value_type*,  host_exec_space> value_type_array_host;

      typedef ExecSpace device_exec_space;
      typedef Kokkos::View<ordinal_type*,device_exec_space> ordinal_type_array_device;
      typedef Kokkos::View<size_type*,   device_exec_space> size_type_array_device;
      typedef Kokkos::View<value_type*,  device_exec_space> value_type_array_device;

      typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

      ///
      /// get the panel dimension of sid
      ///
      inline
      static void
      getSuperPanelSize(const ordinal_type sid,
                        const ordinal_type_array_host &supernodes,
                        const size_type_array_host &sid_super_panel_ptr,
                        const ordinal_type_array_host &blk_super_panel_colidx,
                        /* */ ordinal_type &m,
                        /* */ ordinal_type &n) {
        m = supernodes(sid+1) - supernodes(sid);
        n = blk_super_panel_colidx(sid_super_panel_ptr(sid+1)-1);
      }

      ///
      /// allocate factor space
      ///
      inline
      static void
      allocateSuperPanels(const ordinal_type nsupernodes,
                          const ordinal_type_array_host &supernodes,
                          const size_type_array_host &sid_super_panel_ptr,
                          const ordinal_type_array_host &blk_super_panel_colidx,
                          /* */ size_type_array_host &spanel_ptr,
                          /* */ value_type_array_host &spanel_buf,
                          const ordinal_type_array_host &work) {
        ordinal_type m, n;
        for (ordinal_type sid=0;sid<nsupernodes;++sid) {
          getSuperPanelSize(sid, supernodes, sid_super_panel_ptr, blk_super_panel_colidx, m, n);
          work(sid) = m*n;
        }

        // prefix scan
        spanel_ptr = size_type_array_host("spanel_ptr", nsupernodes+1);
        for (ordinal_type sid=0;sid<nsupernodes;++sid)
          spanel_ptr(sid+1) = spanel_ptr(sid) + work(sid);
        spanel_buf = value_type_array_host("spanel_buf", spanel_ptr(nsupernodes));
      }

      ///
      /// copy sparse matrix to super panels
      ///
      inline
      static void
      copySparseToSuperPanels(// input from sparse matrix
                              const size_type_array_host &ap,
                              const ordinal_type_array_host &aj,
                              const value_type_array_host &ax,
                              const ordinal_type_array_host &perm,
                              const ordinal_type_array_host &peri,
                              // supernodes
                              const ordinal_type nsupernodes,
                              const ordinal_type_array_host &supernodes,
                              const size_type_array_host &gid_super_panel_ptr,
                              const ordinal_type_array_host &gid_super_panel_colidx,
                              const size_type_array_host &sid_super_panel_ptr,
                              const ordinal_type_array_host &blk_super_panel_colidx,
                              // super panel data array
                              const size_type_array_host &super_panel_ptr,
                              const value_type_array_host &super_panel_buf,
                              // work array to store map
                              const ordinal_type_array_host &work) {
        ordinal_type m, n;

        Kokkos::deep_copy(work, -1);        
        
        for (ordinal_type sid=0;sid<nsupernodes;++sid) {
          // get panel for this sid (column major order to use blas and lapack)
          getSuperPanelSize(sid, supernodes, sid_super_panel_ptr, blk_super_panel_colidx, m, n);
          Kokkos::View<value_type**,Kokkos::LayoutLeft,
            typename value_type_array_host::execution_space,
            Kokkos::MemoryUnmanaged> tgt(&super_panel_buf(super_panel_ptr(sid)), m, n);          

          // local to global map
          const ordinal_type goffset = gid_super_panel_ptr(sid);
          for (ordinal_type j=0;j<n;++j) {
            const ordinal_type col = perm(gid_super_panel_colidx(j+goffset)); 
            work(col) = j;
          }

          // row major access to sparse src
          const ordinal_type soffset = supernodes(sid);
          for (ordinal_type i=0;i<m;++i) {
            const ordinal_type row = perm(i+soffset); // row in sparse matrix
            for (ordinal_type k=ap(row);k<ap(row+1);++k) {
              const ordinal_type col = aj(k);
              const ordinal_type j = work(col);
              if (j != -1 && i <= j)  // upper triangular 
                tgt(i, work(col)) = ax(k);
            }
          }
          
          // reset workspace
          for (ordinal_type j=0;j<n;++j) {
            const ordinal_type col = perm(gid_super_panel_colidx(j+goffset)); 
            work(col) = -1;
          }
        }

        // // check copy
        // for (ordinal_type sid=0;sid<nsupernodes;++sid) {
        //   // get panel for this sid (column major order to use blas and lapack)
        //   getSuperPanelSize(sid, supernodes, sid_super_panel_ptr, blk_super_panel_colidx, m, n);
        //   Kokkos::View<value_type**,Kokkos::LayoutLeft,
        //     typename value_type_array_host::execution_space,
        //     Kokkos::MemoryUnmanaged> tgt(&super_panel_buf(super_panel_ptr(sid)), m, n);          

        //   const ordinal_type soffset = supernodes(sid);
        //   const ordinal_type goffset = gid_super_panel_ptr(sid);
        //   for (ordinal_type i=0;i<m;++i) {
        //     for (ordinal_type j=i;j<n;++j) { 
        //       printf("%d %d %f\n", 
        //              (i+soffset), 
        //              (gid_super_panel_colidx(j+goffset)),
        //              tgt(i,j));
        //     }
        //   }
        // }
        // printf("\n\n\n");
      }

    private:
      // matrix input
      ordinal_type _m;
      size_type_array_host _ap;
      ordinal_type_array_host _aj;
      value_type_array_host _ax;

      // graph ordering input
      ordinal_type_array_host _perm, _peri;

      // supernodes input
      ordinal_type _nsupernodes;
      ordinal_type_array_host _supernodes;

      // dof mapping to sparse matrix
      size_type_array_host _gid_super_panel_ptr;
      ordinal_type_array_host _gid_super_panel_colidx;

      // supernode map and panel size configuration
      size_type_array_host _sid_super_panel_ptr;
      ordinal_type_array_host _sid_super_panel_colidx, _blk_super_panel_colidx;

      // supernode tree
      size_type_array_host _stree_ptr;
      ordinal_type_array_host _stree_children, _stree_roots;

      // output : factors
      size_type_array_host _super_panel_ptr;
      value_type_array_host _super_panel_buf;

      ordinal_type_array_host _work;

    private:

      template<typename ViewType>
      KOKKOS_INLINE_FUNCTION
      void
      update(const ordinal_type sid, ViewType &ABR) {
        const size_type sbeg = _sid_super_panel_ptr(sid)+1, send = _sid_super_panel_ptr(sid+1)-1;
        const ordinal_type 
          src_col_beg = _blk_super_panel_colidx(sbeg),
          src_col_end = _blk_super_panel_colidx(send), 
          src_col_size = src_col_end - src_col_beg;

        ordinal_type_array_host map("map", src_col_size);
        const ordinal_type smapoff = _gid_super_panel_ptr(sid);
        auto src_map = Kokkos::subview(_gid_super_panel_colidx, 
                                       range_type(smapoff + src_col_beg,smapoff + src_col_end));
        
        // walk through source rows
        const ordinal_type src_row_offset = _blk_super_panel_colidx(sbeg);
        for (size_type i=sbeg;i<send;++i) {        
          /// ** soruce rows
          const ordinal_type 
            src_row_beg = _blk_super_panel_colidx(i),
            src_row_end = _blk_super_panel_colidx(i+1);

          /// ** target rows
          const ordinal_type row = _sid_super_panel_colidx(i);          
          ordinal_type m, n;
          getSuperPanelSize(row, _supernodes, _sid_super_panel_ptr, _blk_super_panel_colidx, m, n);
          Kokkos::View<value_type**,Kokkos::LayoutLeft,
            typename value_type_array_host::execution_space,
            Kokkos::MemoryUnmanaged> A(&_super_panel_buf(_super_panel_ptr(row)), m, n);          

          /// ** map
          const size_type rbeg = _sid_super_panel_ptr(row), rend = _sid_super_panel_ptr(row+1)-1;
          const ordinal_type 
            tgt_col_beg = _blk_super_panel_colidx(rbeg),
            tgt_col_end = _blk_super_panel_colidx(rend),
            tgt_col_size = tgt_col_end - tgt_col_beg;
          
          const ordinal_type tmapoff = _gid_super_panel_ptr(row);
          auto tgt_map = Kokkos::subview(_gid_super_panel_colidx, 
                                         range_type(tmapoff + tgt_col_beg, tmapoff + tgt_col_end));

          for (ordinal_type k=0,l=0;k<src_col_size;++k) {
            map(k) = -1;
            for (;l<tgt_col_size && tgt_map(l) <= src_map(k);++l)
              if (src_map(k) == tgt_map(l)) { 
                map(k) = l;
                break;
              }
          }
          
          ordinal_type mbeg = 0; for (;map(mbeg) == -1; ++mbeg) ;
          for (ordinal_type jj=mbeg;jj<src_col_size;++jj) {
            const ordinal_type mj = map(jj);
            for (ordinal_type ii=src_row_beg;ii<src_row_end;++ii) {
              const ordinal_type mi = ii+map(mbeg)-src_row_beg;
              A(mi, mj) += ABR(ii-src_row_offset,jj);
            }
          }
        }
      }

      KOKKOS_INLINE_FUNCTION
      void
      recursiveSerialChol(const ordinal_type sid, const ordinal_type sidpar) {
        // recursion
        const ordinal_type ibeg = _stree_ptr(sid), iend = _stree_ptr(sid+1);
        for (ordinal_type i=ibeg;i<iend;++i)
          recursiveSerialChol(_stree_children(i), sid);
        
        ///
        /// body ( panel factorization and herk update
        ///

        // dummy policy and member
        const ordinal_type policy = 0, member = 0;
        
        // chol
        ordinal_type mm, nn;
        getSuperPanelSize(sid, _supernodes, _sid_super_panel_ptr, _blk_super_panel_colidx, mm, nn);
        value_type *ptr = &_super_panel_buf(_super_panel_ptr(sid));

        const ordinal_type m = mm, n = nn - mm;

        Kokkos::View<value_type**,Kokkos::LayoutLeft,
          typename value_type_array_host::execution_space,
          Kokkos::MemoryUnmanaged> ATL(ptr, m, m), ATR(ptr+m*m, m, n);
        
        Chol<Uplo::Upper,Algo::External>
          ::invoke(policy, member, ATL); 

        Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::External>
          ::invoke(policy, member, Diag::Unit(), 1.0, ATL, ATR);

        // if this is not root, it needs to update its parent
        if (sidpar != -1) {
          // temporary workspace ; replaced with memory pool
          Kokkos::View<value_type**,Kokkos::LayoutLeft,
            typename value_type_array_host::execution_space> ABR("ABR", n, n);
          
          Herk<Uplo::Upper,Trans::ConjTranspose,Algo::External>
            ::invoke(policy, member, -1.0, ATR, 1.0, ABR);

          // copy back to its parent
          update(sid, ABR);          
        }
      }

    public:
      NumericTools() = default;
      NumericTools(const NumericTools &b) = default;

      ///
      /// construction (assume input matrix and symbolic are from host)
      ///
      NumericTools(// input matrix A
                   const ordinal_type m,
                   const size_type_array_host &ap,
                   const ordinal_type_array_host &aj,
                   const value_type_array_host &ax,
                   // input permutation
                   const ordinal_type_array_host &perm,
                   const ordinal_type_array_host &peri,
                   // supernodes
                   const ordinal_type nsupernodes,
                   const ordinal_type_array_host &supernodes,
                   const size_type_array_host &gid_super_panel_ptr,
                   const ordinal_type_array_host &gid_super_panel_colidx,
                   const size_type_array_host &sid_super_panel_ptr,
                   const ordinal_type_array_host &sid_super_panel_colidx,
                   const ordinal_type_array_host &blk_super_panel_colidx,
                   const size_type_array_host &stree_ptr,
                   const ordinal_type_array_host &stree_children,
                   const ordinal_type_array_host &stree_roots)
        : _m(m),
          _ap(ap),
          _aj(aj),
          _ax(ax),
          _perm(perm),
          _peri(peri),
          _nsupernodes(nsupernodes),
          _supernodes(supernodes),
          _gid_super_panel_ptr(gid_super_panel_ptr),
          _gid_super_panel_colidx(gid_super_panel_colidx),
          _sid_super_panel_ptr(sid_super_panel_ptr),
          _sid_super_panel_colidx(sid_super_panel_colidx),
          _blk_super_panel_colidx(blk_super_panel_colidx),
          _stree_ptr(stree_ptr),
          _stree_children(stree_children),
          _stree_roots(stree_roots) {}

      ///
      /// host only input (value array can be rewritten in the same sparse structure)
      ///
      inline
      void
      CholeskyFactorize() {
        _work = ordinal_type_array_host("work", _m+1);

        /// allocate for supernode panels
        allocateSuperPanels(_nsupernodes, _supernodes, _sid_super_panel_ptr, _blk_super_panel_colidx,
                            _super_panel_ptr, _super_panel_buf, _work);

        /// copy the input matrix into spanel_buf;
        copySparseToSuperPanels(_ap, _aj, _ax, _perm, _peri,
                                _nsupernodes, _supernodes,
                                _gid_super_panel_ptr, _gid_super_panel_colidx,
                                _sid_super_panel_ptr, _blk_super_panel_colidx,
                                _super_panel_ptr, _super_panel_buf,
                                _work);
        
        const ordinal_type numRoots = _stree_roots.dimension_0();
        for (ordinal_type i=0;i<numRoots;++i) 
          recursiveSerialChol(_stree_roots(i), -1);
      }
      
      // matrix values are only changed (keep workspace)
      inline
      void
      CholeskyFactorize(const value_type_array_host &ax) {
        _ax = ax;

        //Kokkos::deep_copy(_work, 0);
        Kokkos::deep_copy(_super_panel_buf, 0);
        copySparseToSuperPanels(_ap, _aj, _ax, _perm, _peri,
                                _nsupernodes, _supernodes,
                                _gid_super_panel_ptr, _gid_super_panel_colidx,
                                _sid_super_panel_ptr, _blk_super_panel_colidx,
                                _super_panel_ptr, _super_panel_buf,
                                _work);
        

      }
    };

  }
}
#endif
