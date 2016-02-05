#ifndef _KOKKOS_GAUSSSEIDEL_HPP
#define _KOKKOS_GAUSSSEIDEL_HPP

#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_GaussSeidel_impl.hpp"

namespace KokkosKernels{

namespace Experimental{

namespace Graph{

  template <typename KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_>
  void gauss_seidel_symbolic(
      KernelHandle *handle,
      typename KernelHandle::row_lno_t num_rows,
      typename KernelHandle::row_lno_t num_cols,
      lno_row_view_t_ row_map,
      lno_nnz_view_t_ entries,
      bool is_graph_symmetric = true){
    typedef typename Impl::GaussSeidel<KernelHandle, lno_row_view_t_,
          lno_nnz_view_t_, typename KernelHandle::in_scalar_nnz_view_t> SGS;
    SGS sgs(handle,num_rows, num_cols, row_map, entries, is_graph_symmetric);
    sgs.initialize_symbolic();
  }

  template <typename KernelHandle,
            typename lno_row_view_t_,
            typename lno_nnz_view_t_,
            typename scalar_nnz_view_t_>
  void gauss_seidel_numeric(KernelHandle *handle,
      typename KernelHandle::row_lno_t num_rows,
      typename KernelHandle::row_lno_t num_cols,
      lno_row_view_t_ row_map,
      lno_nnz_view_t_ entries,
      scalar_nnz_view_t_ values,
      bool is_graph_symmetric = true
      ){
    typedef typename Impl::GaussSeidel
        <KernelHandle,lno_row_view_t_,
        lno_nnz_view_t_,scalar_nnz_view_t_> SGS;
    SGS sgs(handle, num_rows, num_cols, row_map, entries, values, is_graph_symmetric);
    sgs.initialize_numeric();
  }

  template <typename KernelHandle,
    typename lno_row_view_t_,
    typename lno_nnz_view_t_,
    typename scalar_nnz_view_t_,
    typename x_scalar_view_t,
    typename y_scalar_view_t>
  void symmetric_gauss_seidel_apply(KernelHandle *handle,
      typename KernelHandle::row_lno_t num_rows,
      typename KernelHandle::row_lno_t num_cols,
      lno_row_view_t_ row_map,
      lno_nnz_view_t_ entries,
      scalar_nnz_view_t_ values,
      x_scalar_view_t x_lhs_output_vec,
      y_scalar_view_t y_rhs_input_vec,
      bool init_zero_x_vector = false,
      bool update_y_vector = true,
      int numIter = 1){
    typedef
        typename Impl::GaussSeidel
            <KernelHandle,
            lno_row_view_t_, lno_nnz_view_t_,scalar_nnz_view_t_ > SGS;

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
  typename lno_row_view_t_,
  typename lno_nnz_view_t_,
  typename scalar_nnz_view_t_,
  typename x_scalar_view_t, typename y_scalar_view_t>
  void forward_sweep_gauss_seidel_apply(
      KernelHandle *handle,
      typename KernelHandle::row_lno_t num_rows,
      typename KernelHandle::row_lno_t num_cols,
      lno_row_view_t_ row_map,
      lno_nnz_view_t_ entries,
      scalar_nnz_view_t_ values,
      x_scalar_view_t x_lhs_output_vec,
      y_scalar_view_t y_rhs_input_vec,
      bool init_zero_x_vector = false,
      bool update_y_vector = true,
      int numIter = 1){
    typedef typename Impl::GaussSeidel <KernelHandle,
            lno_row_view_t_, lno_nnz_view_t_,scalar_nnz_view_t_ > SGS;

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
  typename lno_row_view_t_,
  typename lno_nnz_view_t_,
  typename scalar_nnz_view_t_,
  typename x_scalar_view_t, typename y_scalar_view_t>
  void backward_sweep_gauss_seidel_apply(
      KernelHandle *handle,
      typename KernelHandle::row_lno_t num_rows,
      typename KernelHandle::row_lno_t num_cols,
      lno_row_view_t_ row_map,
      lno_nnz_view_t_ entries,
      scalar_nnz_view_t_ values,
      x_scalar_view_t x_lhs_output_vec,
      y_scalar_view_t y_rhs_input_vec,
      bool init_zero_x_vector = false,
      bool update_y_vector = true,
      int numIter = 1){
    typedef typename Impl::GaussSeidel <KernelHandle,
            lno_row_view_t_, lno_nnz_view_t_,scalar_nnz_view_t_ > SGS;
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
