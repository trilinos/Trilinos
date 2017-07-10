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

      typedef typename crs_matrix_type::ordinal_type_array ordinal_type_array;
      typedef typename crs_matrix_type::size_type_array size_type_array;
      typedef typename crs_matrix_type::value_type_array value_type_array;

      typedef Kokkos::View<value_type**,Kokkos::LayoutLeft,exec_space> value_type_matrix;

      typedef Kokkos::Future<int,exec_space> future_type;
      typedef Kokkos::View<future_type*,exec_space> future_type_array;

      ///
      /// Phase 1: symbolic
      ///

      // supernodes input
      ConstUnmanagedViewType<ordinal_type_array> supernodes;

      // dof mapping to sparse matrix
      ConstUnmanagedViewType<size_type_array> gid_super_panel_ptr;
      ConstUnmanagedViewType<ordinal_type_array> gid_super_panel_colidx;

      // supernode map and panel size configuration
      ConstUnmanagedViewType<size_type_array> sid_super_panel_ptr;
      ConstUnmanagedViewType<ordinal_type_array> sid_super_panel_colidx, blk_super_panel_colidx;

      // supernode tree
      ConstUnmanagedViewType<size_type_array> stree_ptr;
      ConstUnmanagedViewType<ordinal_type_array> stree_children;

      ///
      /// Phase 2: factors
      ///
      ConstUnmanagedViewType<size_type_array> super_panel_ptr;
      UnmanagedViewType<value_type_array> super_panel_buf;

      // ///
      // /// Phase 3: abr (schur complements)
      // ///
      ordinal_type max_supernode_size, max_schur_size;

      ///
      /// Phase 4: solve (rhs multivector)
      UnmanagedViewType<value_type_matrix> x;

      KOKKOS_INLINE_FUNCTION
      SupernodeInfo() = default;

      KOKKOS_INLINE_FUNCTION
      SupernodeInfo(const SupernodeInfo &b) = default;

      KOKKOS_INLINE_FUNCTION
      void
      getSuperPanelSize(const ordinal_type sid,
                        /* */ ordinal_type &m,
                        /* */ ordinal_type &n) const {
        m = supernodes(sid+1) - supernodes(sid);
        n = blk_super_panel_colidx(sid_super_panel_ptr(sid+1)-1);
      }

      KOKKOS_INLINE_FUNCTION
      void
      getSuperPanel(const ordinal_type sid,
                    const ordinal_type m,
                    const ordinal_type n,
                    /* */ UnmanagedViewType<value_type_matrix> &A) const {
        A = value_type_matrix(&super_panel_buf(super_panel_ptr(sid)), m, n);
      }

      KOKKOS_INLINE_FUNCTION
      value_type*
      getSuperPanelPtr(const ordinal_type sid) const {
        return &super_panel_buf(super_panel_ptr(sid));
      }

      inline
      void
      allocateSuperPanels(/* */ size_type_array &spanel_ptr,
                          /* */ value_type_array &spanel_buf,
                          const ordinal_type_array &work) {

        max_schur_size = 0; 
        max_supernode_size = 0;
        const ordinal_type nsupernodes = supernodes.dimension_0() - 1;
        for (ordinal_type sid=0;sid<nsupernodes;++sid) {
          ordinal_type m, n;
          getSuperPanelSize(sid, m, n);
          work(sid) = m*n;
          max_supernode_size = max(max_supernode_size, m);
          max_schur_size = max(max_schur_size, n-m);
        }

        // prefix scan
        spanel_ptr = size_type_array("super_panel_ptr", nsupernodes+1);
        for (ordinal_type sid=0;sid<nsupernodes;++sid)
          spanel_ptr(sid+1) = spanel_ptr(sid) + work(sid);
        spanel_buf = value_type_array("super_panel_buf", spanel_ptr(nsupernodes));
      }

      // inline
      // void
      // allocateWorkspacePerSupernode(/* */ size_type_array &schur_ptr,
      //                               /* */ value_type_array &schur_buf,
      //                               const ordinal_type_array &work) {
      //   const ordinal_type nsupernodes = supernodes.dimension_0() - 1;
      //   for (ordinal_type sid=0;sid<nsupernodes;++sid) {
      //     // workspace for schur complement
      //     ordinal_type m, n;
      //     getSuperPanelSize(sid, m, n);

      //     // workspace for update map
      //     const ordinal_type
      //       cbeg = blk_super_panel_colidx(sid_super_panel_ptr(sid)+1),
      //       cend = blk_super_panel_colidx(sid_super_panel_ptr(sid+1)-1);

      //     work(sid) = (m-n)*(m-n) + (cend - cbeg);
      //   }

      //   // prefix scan
      //   schur_ptr = size_type_array("super_schur_ptr", nsupernodes+1);
      //   for (ordinal_type sid=0;sid<nsupernodes;++sid)
      //     schur_ptr(sid+1) = schur_ptr(sid) + work(sid);
      //   schur_buf = value_type_array("super_schur_buf", schur_ptr(nsupernodes) + 1);
      // }

      // inline
      // ordinal_type
      // computeMaxSchurDimension() {
      //   const ordinal_type nsupernodes = supernodes.dimension_0() - 1;
      //   ordinal_type max_schur_dimension = 0;
      //   for (ordinal_type sid=0;sid<nsupernodes;++sid) {
      //     ordinal_type m, n;
      //     getSuperPanelSize(sid, m, n);
      //     max_schur_dimension = max(max_schur_dimension, (n-m));
      //   }
      //   return max_schur_dimension;
      // }

      // inline
      // size_type
      // computeWorkspaceSerialChol() {
      //   const ordinal_type nsupernodes = supernodes.dimension_0() - 1;
      //   size_type workspace = 0;
      //   for (ordinal_type sid=0;sid<nsupernodes;++sid) {
      //     ordinal_type m, n;
      //     getSuperPanelSize(sid, m, n);
      //     workspace = max(workspace, (n-m)*(n-m));
      //   }
      //   return workspace;
      // }

      // inline
      // size_type
      // computeWorkspaceCholByBlocks(const ordinal_type mb) {
      //   const ordinal_type nsupernodes = supernodes.dimension_0() - 1;
      //   size_type workspace = 0;
      //   for (ordinal_type sid=0;sid<nsupernodes;++sid) {
      //     ordinal_type m, n;
      //     getSuperPanelSize(sid, m, n);

      //     const ordinal_type bm = m/mb + 1, bn = (n-m)/mb + 1;
      //     workspace = max(workspace, bm*bm + bm*bn + bn*bn);
      //   }
      //   return workspace;
      // }

      inline
      void
      copySparseToSuperPanels(// input from sparse matrix
                              const size_type_array &ap,
                              const ordinal_type_array &aj,
                              const value_type_array &ax,
                              const ordinal_type_array &perm,
                              const ordinal_type_array &peri,
                              // work array to store map
                              const ordinal_type_array &work) {
        const ordinal_type nsupernodes = supernodes.dimension_0() - 1;

        Kokkos::deep_copy(work, -1);

        for (ordinal_type sid=0;sid<nsupernodes;++sid) {
          // grab super panel
          ordinal_type m, n;
          getSuperPanelSize(sid, m, n);

          UnmanagedViewType<value_type_matrix>
            tgt(&super_panel_buf(super_panel_ptr(sid)), m, n);

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
      }

      inline
      crs_matrix_type
      createCrsMatrix(const bool replace_value_with_one = false) {
        // count m, n, nnz
        const ordinal_type
          nsupernodes = supernodes.dimension_0() - 1,
          mm = supernodes(nsupernodes),
          nn = mm;

        size_type cnt = 0;
        size_type_array ap("ap", mm+1);
        for (ordinal_type sid=0;sid<nsupernodes;++sid) {
          // grab super panel
          ordinal_type m, n;
          getSuperPanelSize(sid, m, n);

          // row major access to sparse src
          const ordinal_type soffset = supernodes(sid);
          for (ordinal_type i=0;i<m;++i) {
            const ordinal_type row = i+soffset; // row in sparse matrix
            ap(row) = cnt;
            cnt += (n - i); // upper triangular only
          }
        }
        ap(mm) = cnt;

        // fill the matrix
        const size_type nnz = cnt; cnt = 0;
        ordinal_type_array aj("aj", nnz);
        value_type_array ax("ax", nnz);
        for (ordinal_type sid=0;sid<nsupernodes;++sid) {

          // grab super panel
          ordinal_type m, n;
          getSuperPanelSize(sid, m, n);

          UnmanagedViewType<value_type_matrix>
            src(&super_panel_buf(super_panel_ptr(sid)), m, n);

          // row major access to sparse src
          const ordinal_type
            soffset = supernodes(sid),
            goffset = gid_super_panel_ptr(sid);

          for (ordinal_type i=0;i<m;++i) {
            const size_type beg = ap(i+soffset);
            for (ordinal_type j=i,k=0;j<n;++j,++k) {
              const ordinal_type col = gid_super_panel_colidx(j+goffset);
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
