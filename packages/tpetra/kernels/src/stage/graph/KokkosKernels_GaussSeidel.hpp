#ifndef _KOKKOS_GAUSSSEIDEL_HPP
#define _KOKKOS_GAUSSSEIDEL_HPP

#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_GaussSeidel_impl.hpp"

namespace KokkosKernels{

namespace Experimental{

namespace Graph{

  template <typename KernelHandle, typename in_row_index_view_type, typename in_nonzero_index_view_type>
  void gauss_seidel_symbolic(
      KernelHandle *handle,
      typename KernelHandle::row_index_type num_rows,
      typename KernelHandle::row_index_type num_cols,
      in_row_index_view_type row_map,
      in_nonzero_index_view_type entries,
      bool is_graph_symmetric = true){
    typedef typename Impl::GaussSeidel<KernelHandle, in_row_index_view_type,
          in_nonzero_index_view_type, typename KernelHandle::in_nonzero_value_view_type> SGS;
    SGS sgs(handle,num_rows, num_cols, row_map, entries, is_graph_symmetric);
    sgs.initialize_symbolic();
  }

  template <typename KernelHandle,
            typename in_row_index_view_type,
            typename in_nonzero_index_view_type,
            typename in_nonzero_value_view_type>
  void gauss_seidel_numeric(KernelHandle *handle,
      typename KernelHandle::row_index_type num_rows,
      typename KernelHandle::row_index_type num_cols,
      in_row_index_view_type row_map,
      in_nonzero_index_view_type entries,
      in_nonzero_value_view_type values,
      bool is_graph_symmetric = true
      ){
    typedef typename Impl::GaussSeidel
        <KernelHandle,in_row_index_view_type,
        in_nonzero_index_view_type,in_nonzero_value_view_type> SGS;
    SGS sgs(handle, num_rows, num_cols, row_map, entries, values, is_graph_symmetric);
    sgs.initialize_numeric();
  }

  template <typename KernelHandle,
    typename in_row_index_view_type,
    typename in_nonzero_index_view_type,
    typename in_nonzero_value_view_type,
    typename x_value_array_type,
    typename y_value_array_type>
  void symmetric_gauss_seidel_apply(KernelHandle *handle,
      typename KernelHandle::row_index_type num_rows,
      typename KernelHandle::row_index_type num_cols,
      in_row_index_view_type row_map,
      in_nonzero_index_view_type entries,
      in_nonzero_value_view_type values,
      x_value_array_type x_lhs_output_vec,
      y_value_array_type y_rhs_input_vec,
      bool init_zero_x_vector = false,
      bool update_y_vector = true,
      int numIter = 1){
    typedef
        typename Impl::GaussSeidel
            <KernelHandle,
            in_row_index_view_type, in_nonzero_index_view_type,in_nonzero_value_view_type > SGS;

    SGS sgs(handle, num_rows, num_cols, row_map, entries, values);
    sgs.apply(
        x_lhs_output_vec,
        y_rhs_input_vec,
        init_zero_x_vector,
        numIter,
        true,
        true,update_y_vector);

  }

  template <class KernelHandle,
  typename in_row_index_view_type,
  typename in_nonzero_index_view_type,
  typename in_nonzero_value_view_type,
  typename x_value_array_type, typename y_value_array_type>
  void forward_sweep_gauss_seidel_apply(
      KernelHandle *handle,
      typename KernelHandle::row_index_type num_rows,
      typename KernelHandle::row_index_type num_cols,
      in_row_index_view_type row_map,
      in_nonzero_index_view_type entries,
      in_nonzero_value_view_type values,
      x_value_array_type x_lhs_output_vec,
      y_value_array_type y_rhs_input_vec,
      bool init_zero_x_vector = false,
      bool update_y_vector = true,
      int numIter = 1){
    typedef typename Impl::GaussSeidel <KernelHandle,
            in_row_index_view_type, in_nonzero_index_view_type,in_nonzero_value_view_type > SGS;

    SGS sgs(handle, num_rows, num_cols, row_map, entries, values);
    sgs.apply(
        x_lhs_output_vec,
        y_rhs_input_vec,
        init_zero_x_vector,
        numIter,
        true,
        false,update_y_vector);

  }
  template <class KernelHandle,
  typename in_row_index_view_type,
  typename in_nonzero_index_view_type,
  typename in_nonzero_value_view_type,
  typename x_value_array_type, typename y_value_array_type>
  void backward_sweep_gauss_seidel_apply(
      KernelHandle *handle,
      typename KernelHandle::row_index_type num_rows,
      typename KernelHandle::row_index_type num_cols,
      in_row_index_view_type row_map,
      in_nonzero_index_view_type entries,
      in_nonzero_value_view_type values,
      x_value_array_type x_lhs_output_vec,
      y_value_array_type y_rhs_input_vec,
      bool init_zero_x_vector = false,
      bool update_y_vector = true,
      int numIter = 1){
    typedef typename Impl::GaussSeidel <KernelHandle,
            in_row_index_view_type, in_nonzero_index_view_type,in_nonzero_value_view_type > SGS;
    SGS sgs(handle, num_rows, num_cols, row_map, entries, values);
    sgs.apply(
        x_lhs_output_vec,
        y_rhs_input_vec,
        init_zero_x_vector,
        numIter,
        false,
        true, update_y_vector);

  }
}
}
}
#endif
