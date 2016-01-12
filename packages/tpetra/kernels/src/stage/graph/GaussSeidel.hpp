#ifndef _KOKKOS_GAUSSSEIDEL_HPP
#define _KOKKOS_GAUSSSEIDEL_HPP

#include "KokkosKernelsHandle.hpp"
#include "GaussSeidel_impl.hpp"

namespace KokkosKernels{

namespace Experimental{

namespace Graph{

  template <class KernelHandle>
  void gauss_seidel_symbolic(
      KernelHandle *handle,
      typename KernelHandle::idx num_rows,
      typename KernelHandle::idx num_cols,
      typename KernelHandle::idx_array_type row_map,
      typename KernelHandle::idx_edge_array_type entries,
      bool is_graph_symmetric = true){
    typedef typename Impl::GaussSeidel<KernelHandle> SGS;
    SGS sgs(handle,num_rows, num_cols, row_map, entries, is_graph_symmetric);
    sgs.initialize_symbolic();
  }

  template <class KernelHandle>
  void gauss_seidel_numeric(KernelHandle *handle,
      typename KernelHandle::idx num_rows,
      typename KernelHandle::idx num_cols,
      typename KernelHandle::idx_array_type row_map,
      typename KernelHandle::idx_edge_array_type entries,
      typename KernelHandle::value_array_type values,
      bool is_graph_symmetric = true
      ){
    typedef typename Impl::GaussSeidel<KernelHandle> SGS;
    SGS sgs(handle, num_rows, num_cols, row_map, entries, values, is_graph_symmetric);
    sgs.initialize_numeric();
  }

  template <class KernelHandle>
  void symmetric_gauss_seidel_apply(KernelHandle *handle,
      typename KernelHandle::idx num_rows,
      typename KernelHandle::idx num_cols,
      typename KernelHandle::idx_array_type row_map,
      typename KernelHandle::idx_edge_array_type entries,
      typename KernelHandle::value_array_type values,
      typename KernelHandle::value_array_type x_lhs_output_vec,
      typename KernelHandle::value_array_type y_rhs_input_vec,
      bool init_zero_x_vector = false,
      int numIter = 1){
    typedef typename Impl::GaussSeidel<KernelHandle> SGS;
    SGS sgs(handle, num_rows, num_cols, row_map, entries, values);
    sgs.apply(
        x_lhs_output_vec,
        y_rhs_input_vec,
        init_zero_x_vector,
        numIter,
        true,
        true);

  }

  template <class KernelHandle>
  void forward_sweep_gauss_seidel_apply(KernelHandle *handle,
      typename KernelHandle::idx num_rows,
      typename KernelHandle::idx num_cols,
      typename KernelHandle::idx_array_type row_map,
      typename KernelHandle::idx_edge_array_type entries,
      typename KernelHandle::value_array_type values,
      typename KernelHandle::value_array_type x_lhs_output_vec,
      typename KernelHandle::value_array_type y_rhs_input_vec,
      bool init_zero_x_vector = false,
      int numIter = 1){
    typedef typename Impl::GaussSeidel<KernelHandle> SGS;
    SGS sgs(handle, num_rows, num_cols, row_map, entries, values);
    sgs.apply(
        x_lhs_output_vec,
        y_rhs_input_vec,
        init_zero_x_vector,
        numIter,
        true,
        false);

  }
  template <class KernelHandle>
  void backward_sweep_gauss_seidel_apply(KernelHandle *handle,
      typename KernelHandle::idx num_rows,
      typename KernelHandle::idx num_cols,
      typename KernelHandle::idx_array_type row_map,
      typename KernelHandle::idx_edge_array_type entries,
      typename KernelHandle::value_array_type values,
      typename KernelHandle::value_array_type x_lhs_output_vec,
      typename KernelHandle::value_array_type y_rhs_input_vec,
      bool init_zero_x_vector = false,
      int numIter = 1){
    typedef typename Impl::GaussSeidel<KernelHandle> SGS;
    SGS sgs(handle, num_rows, num_cols, row_map, entries, values);
    sgs.apply(
        x_lhs_output_vec,
        y_rhs_input_vec,
        init_zero_x_vector,
        numIter,
        false,
        true);

  }
}
}
}
#endif
