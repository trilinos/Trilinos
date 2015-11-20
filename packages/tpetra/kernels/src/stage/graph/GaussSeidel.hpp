#ifndef _KOKKOS_GAUSSSEIDEL_HPP
#define _KOKKOS_GAUSSSEIDEL_HPP

#include "KokkosKernelsHandle.hpp"
#include "GaussSeidel_impl.hpp"
namespace Experimental{

namespace KokkosKernels{
namespace Graph{

  template <class KernelHandle>
  void gauss_seidel_symbolic(KernelHandle *handle,
      typename KernelHandle::idx_array_type row_map,
      typename KernelHandle::idx_edge_array_type entries){
    typedef typename Experimental::KokkosKernels::Graph::Impl::GaussSeidel<KernelHandle> SGS;
    SGS sgs(handle, row_map, entries);
    sgs.initialize_symbolic();
  }

  template <class KernelHandle>
  void gauss_seidel_numeric(KernelHandle *handle,
      typename KernelHandle::idx_array_type row_map,
      typename KernelHandle::idx_edge_array_type entries,
      typename KernelHandle::value_array_type values
      ){
    typedef typename Experimental::KokkosKernels::Graph::Impl::GaussSeidel<KernelHandle> SGS;
    SGS sgs(handle, row_map, entries, values);
    sgs.initialize_numeric();
  }

  template <class KernelHandle>
  void gauss_seidel_apply(KernelHandle *handle,
      typename KernelHandle::value_array_type &x_lhs_output_vec,
      typename KernelHandle::value_array_type y_rhs_input_vec,
      bool init_zero_x_vector = false,
      int numIter = 1){
    typedef typename Experimental::KokkosKernels::Graph::Impl::GaussSeidel<KernelHandle> SGS;
    SGS sgs(handle);
    sgs.apply(
        x_lhs_output_vec,
        y_rhs_input_vec,
        init_zero_x_vector,
        numIter);

  }

}
}
}
#endif
