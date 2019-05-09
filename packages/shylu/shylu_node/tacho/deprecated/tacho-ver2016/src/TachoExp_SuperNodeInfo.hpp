#ifndef __TACHOEXP_SUPERNODE_INFO_HPP__
#define __TACHOEXP_SUPERNODE_INFO_HPP__

#include "TachoExp_Util.hpp"

/// \file TachoExp_SuperNodeInfo.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  namespace Experimental {

    template<typename ValueType, typename ExecSpace>
    struct SuperNodeInfo {
      typedef ValueType value_type;
      typedef ExecSpace exec_space;

      typedef Kokkos::View<ordinal_type*,exec_space> ordinal_type_array;
      typedef Kokkos::View<size_type*,   exec_space> size_type_array;
      typedef Kokkos::View<value_type*,  exec_space> value_type_array;

      typedef Kokkos::View<value_type**,Kokkos::LayoutLeft,exec_space> value_type_matrix;

      typedef Kokkos::Future<int,exec_space> future_type;
      typedef Kokkos::View<future_type*,exec_space> future_type_array;      

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

      // output : factors
      ConstUnmanagedViewType<size_type_array> super_panel_ptr;
      UnmanagedViewType<value_type_array> super_panel_buf;


      // parallel version; fugure list corresponding to each supernode 
      // temporal memory allocation depends on kokkos memory pool
      //UnmanagedViewType<future_type_array> supernodes_future;
      future_type_array supernodes_future;

      // static work space for serial execution (not sure if this is necessary)
      UnmanagedViewType<value_type_array> super_panel_serial_work;

      KOKKOS_INLINE_FUNCTION
      SuperNodeInfo() = default;

      KOKKOS_INLINE_FUNCTION
      SuperNodeInfo(const SuperNodeInfo &b) = default;

      KOKKOS_INLINE_FUNCTION
      void
      getSuperPanelSize(const ordinal_type sid,
                        /* */ ordinal_type &m,
                        /* */ ordinal_type &n) {
        m = supernodes(sid+1) - supernodes(sid);
        n = blk_super_panel_colidx(sid_super_panel_ptr(sid+1)-1);
      }

      KOKKOS_INLINE_FUNCTION
      void
      getSuperPanel(const ordinal_type sid, 
                    const ordinal_type m, 
                    const ordinal_type n, 
                    /* */ UnmanagedViewType<value_type_matrix> &A) {
        A = value_type_matrix(&super_panel_buf(super_panel_ptr(sid)), m, n);  
      }

      KOKKOS_INLINE_FUNCTION
      value_type*
      getSuperPanelPtr(const ordinal_type sid) {
        return &super_panel_buf(super_panel_ptr(sid));
      }
    };

  }
}

#endif
