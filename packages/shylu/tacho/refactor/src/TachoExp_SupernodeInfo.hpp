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

      typedef Kokkos::View<ordinal_type*,exec_space> ordinal_type_array;
      typedef Kokkos::View<size_type*,   exec_space> size_type_array;
      typedef Kokkos::View<value_type*,  exec_space> value_type_array;

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

      ///
      /// Phase 3: abr (schur complements)
      ///
      ConstUnmanagedViewType<size_type_array> super_schur_ptr;
      UnmanagedViewType<value_type_array> super_schur_buf;
      
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
        const ordinal_type nsupernodes = supernodes.dimension_0() - 1;
        for (ordinal_type sid=0;sid<nsupernodes;++sid) {
          ordinal_type m, n;
          getSuperPanelSize(sid, m, n);
          work(sid) = m*n;
        }

        // prefix scan
        spanel_ptr = size_type_array("super_panel_ptr", nsupernodes+1);
        for (ordinal_type sid=0;sid<nsupernodes;++sid)
          spanel_ptr(sid+1) = spanel_ptr(sid) + work(sid);
        spanel_buf = value_type_array("super_panel_buf", spanel_ptr(nsupernodes));
      }

      inline
      void
      allocateWorkspacePerSupernode(/* */ size_type_array &schur_ptr,
                                    /* */ value_type_array &schur_buf,
                                    const ordinal_type_array &work) {
        const ordinal_type nsupernodes = supernodes.dimension_0() - 1;
        for (ordinal_type sid=0;sid<nsupernodes;++sid) {
          ordinal_type m, n;
          getSuperPanelSize(sid, m, n);
          work(sid) = (m-n)*(m-n);
        }

        // prefix scan
        schur_ptr = size_type_array("super_schur_ptr", nsupernodes+1);
        for (ordinal_type sid=0;sid<nsupernodes;++sid)
          schur_ptr(sid+1) = schur_ptr(sid) + work(sid);
        schur_buf = value_type_array("super_schur_buf", schur_ptr(nsupernodes));
      }

      inline
      size_type
      computeWorkspaceSerialChol() {
        const ordinal_type nsupernodes = supernodes.dimension_0() - 1;
        size_type workspace = 0;
        for (ordinal_type sid=0;sid<nsupernodes;++sid) {
          // supernodes are myself, parent, empty one (range is used for blocks it requires end point)
          const bool is_direct_update = (sid_super_panel_ptr(sid+1) - sid_super_panel_ptr(sid)) == 3;
          if (!is_direct_update) {
            ordinal_type m, n;
            getSuperPanelSize(sid, m, n);
            workspace = max(workspace, (n-m)*(n-m));
          }
        }
        return workspace;
      }

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
    };

  }
}

#endif
