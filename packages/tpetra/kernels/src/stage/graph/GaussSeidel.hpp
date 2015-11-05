#ifndef _KOKKOS_GAUSSSEIDEL_HPP
#define _KOKKOS_GAUSSSEIDEL_HPP

#include "KokkosKernelsHandle.hpp"
#include "GaussSeidel_impl.hpp"
namespace Experimental{

namespace KokkosKernels{
namespace Graph{
  template <class KernelHandle>
  void init_gauss_seidel_symbolic(KernelHandle *handle){
    typedef typename Experimental::KokkosKernels::Graph::Impl::GaussSeidel<KernelHandle> SGS;
    SGS sgs(handle);
    sgs.initialize_symbolic();
  }
  template <class KernelHandle>
  void init_gauss_seidel_numeric(KernelHandle *handle){
    typedef typename Experimental::KokkosKernels::Graph::Impl::GaussSeidel<KernelHandle> SGS;
    SGS sgs(handle);
    sgs.initialize_numeric();
  }
  template <class KernelHandle>
  void init_gauss_seidel_solve(KernelHandle *handle){
    typedef typename Experimental::KokkosKernels::Graph::Impl::GaussSeidel<KernelHandle> SGS;
    SGS sgs(handle);
    sgs.initialize_solve();
  }

  template <class KernelHandle>
  void apply_gauss_seidel_solve(KernelHandle *handle,
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

  template <class KernelHandle>
  void apply_gauss_seidel_numeric(KernelHandle *handle,
      typename KernelHandle::value_array_type &x_lhs_output_vec,
      typename KernelHandle::value_array_type y_rhs_input_vec,
      bool init_zero_x_vector = false,
      int numIter = 1){
    apply_gauss_seidel_solve(handle, x_lhs_output_vec, y_rhs_input_vec, init_zero_x_vector, numIter);

  }

  template <class KernelHandle>
  void apply_gauss_seidel_symbolic(KernelHandle *handle,
      typename KernelHandle::value_array_type &x_lhs_output_vec,
      typename KernelHandle::value_array_type y_rhs_input_vec,
      bool init_zero_x_vector = false,
      int numIter = 1){
    apply_gauss_seidel_solve(handle, x_lhs_output_vec, y_rhs_input_vec, init_zero_x_vector, numIter);

  }

}
}
}
#endif
