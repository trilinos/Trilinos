
#ifndef _KOKKOSSPGEMMCUSP_HPP
#define _KOKKOSSPGEMMCUSP_HPP

#ifdef KERNELS_HAVE_CUSP
#include <cusp/multiply.h>
#include <cusp/csr_matrix.h>
#endif

namespace KokkosKernels{

namespace Experimental{

namespace Graph{
namespace Impl{


template <typename cusparray, typename kokkosarray>
struct CopyArrayToCuspArray{
  cusparray c;
  kokkosarray *k;

  CopyArrayToCuspArray(cusparray &c_, kokkosarray *k_): c(c_), k(k_){}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &i) const {
    c[i] = k[i];
  }
};




template <typename KernelHandle,
  typename in_row_index_view_type,
  typename in_nonzero_index_view_type,
  typename in_nonzero_value_view_type>
void CUSP_apply(
    KernelHandle *handle,
    typename KernelHandle::row_lno_t m,
    typename KernelHandle::row_lno_t n,
    typename KernelHandle::row_lno_t k,
    in_row_index_view_type row_mapA,
    in_nonzero_index_view_type entriesA,
    in_nonzero_value_view_type valuesA,

    bool transposeA,
    in_row_index_view_type row_mapB,
    in_nonzero_index_view_type entriesB,
    in_nonzero_value_view_type valuesB,
    bool transposeB,
    typename in_row_index_view_type::non_const_type &row_mapC,
    typename in_nonzero_index_view_type::non_const_type &entriesC,
    typename in_nonzero_value_view_type::non_const_type &valuesC){
#ifdef KERNELS_HAVE_CUSP
  typedef typename KernelHandle::row_lno_t idx;
  typedef typename KernelHandle::nnz_scalar_t value_type;

  typedef typename in_row_index_view_type::device_type device1;
  typedef typename in_nonzero_index_view_type::device_type device2;
  typedef typename in_nonzero_value_view_type::device_type device3;

  std::cout << "RUNNING CUSP" << std::endl;
  if (Kokkos::Impl::is_same<Kokkos::Cuda, device1 >::value){
    std::cerr << "MEMORY IS NOT ALLOCATED IN GPU DEVICE for CUSP" << std::endl;
    return;
  }
  if (Kokkos::Impl::is_same<Kokkos::Cuda, device2 >::value){
    std::cerr << "MEMORY IS NOT ALLOCATED IN GPU DEVICE for CUSP" << std::endl;
    return;
  }
  if (Kokkos::Impl::is_same<Kokkos::Cuda, device3 >::value){
    std::cerr << "MEMORY IS NOT ALLOCATED IN GPU DEVICE for CUSP" << std::endl;
    return;
  }

  typedef in_row_index_view_type idx_array_type;

  typedef typename Kokkos::RangePolicy<typename KernelHandle::HandleExecSpace> my_exec_space;




  idx nnzA = entriesA.dimension_0();
  idx nnzB = entriesB.dimension_0();


  idx *a_xadj = (idx *)row_mapA.ptr_on_device();
  idx *b_xadj = (idx *)row_mapB.ptr_on_device();


  idx *a_adj = (idx *)entriesA.ptr_on_device();
  idx *b_adj = (idx *)entriesB.ptr_on_device();

  value_type *a_ew = valuesA.ptr_on_device();
  value_type *b_ew = valuesB.ptr_on_device();

  /*
  thrust::device_ptr<idx> dev_a_xadj(a_xadj);
  thrust::device_ptr<idx> dev_a_adj(a_adj);
  thrust::device_ptr<idx> dev_b_xadj(b_xadj);
  thrust::device_ptr<idx> dev_b_adj(b_adj);
  thrust::device_ptr<value_type> dev_a_ew(a_ew);
  thrust::device_ptr<value_type> dev_b_ew(b_ew);

  */
  typedef typename cusp::array1d_view< thrust::device_ptr<idx> > IDXArray1dView;
  typedef typename cusp::array1d_view< thrust::device_ptr<value_type> > VALUEArray1dView;
  //typedef typename cusp::array1d<idx, cusp::device_memory> IDXArray1dView;
  //typedef typename cusp::array1d<value_type, cusp::device_memory> VALUEArray1dView;
  IDXArray1dView arraya_xadj(thrust::device_pointer_cast(a_xadj), thrust::device_pointer_cast(a_xadj) + m + 1);
  IDXArray1dView arraya_adj(thrust::device_pointer_cast(a_adj), thrust::device_pointer_cast(a_adj) + nnzA);
  IDXArray1dView arrayb_xadj(thrust::device_pointer_cast(b_xadj), thrust::device_pointer_cast(b_xadj) + n + 1);
  IDXArray1dView arrayb_adj(thrust::device_pointer_cast(b_adj), thrust::device_pointer_cast(b_adj) + nnzB);
  VALUEArray1dView arraya_ew(thrust::device_pointer_cast(a_ew), thrust::device_pointer_cast(a_ew) + nnzA);
  VALUEArray1dView arrayb_ew(thrust::device_pointer_cast(b_ew), thrust::device_pointer_cast(b_ew)+ nnzB);

  typedef typename cusp::csr_matrix_view<IDXArray1dView, IDXArray1dView, VALUEArray1dView, idx,value_type,cusp::device_memory> cuspMatrix_View;

  cuspMatrix_View A(m, n, entriesA.dimension_0(), arraya_xadj, arraya_adj, arraya_ew);
  cuspMatrix_View B(n, k, entriesB.dimension_0(), arrayb_xadj, arrayb_adj, arrayb_ew);

  /*
  CopyArrayToCuspArray<typename cuspMatrix::row_offsets_array_type, typename KernelHandle::idx_array_type> Aforward(A.row_offsets, row_mapA);
  Kokkos::parallel_for (my_exec_space (0, m + 1) , Aforward);
  Kokkos::parallel_for (my_exec_space (0, n + 1) , CopyArrayToCuspArray<typename cuspMatrix::row_offsets_array_type, typename KernelHandle::idx_array_type>(B.row_offsets, row_mapB));

  Kokkos::parallel_for (my_exec_space (0, entriesA.dimension_0()) , CopyArrayToCuspArray<typename cuspMatrix::column_indices_array_type, typename KernelHandle::idx_edge_array_type>(A.column_indices, entriesA));
  Kokkos::parallel_for (my_exec_space (0, entriesB.dimension_0()) , CopyArrayToCuspArray<typename cuspMatrix::column_indices_array_type, typename KernelHandle::idx_edge_array_type>(B.column_indices, entriesB));

  Kokkos::parallel_for (my_exec_space (0, valuesA.dimension_0()) , CopyArrayToCuspArray<typename cuspMatrix::values_array_type, typename KernelHandle::value_array_type>(A.values, valuesA));
  Kokkos::parallel_for (my_exec_space (0, valuesB.dimension_0()) , CopyArrayToCuspArray<typename cuspMatrix::values_array_type, typename KernelHandle::value_array_type>(B.values, valuesB));
  */

  typedef typename cusp::csr_matrix<idx,value_type,cusp::device_memory> cuspMatrix;
  //typedef cuspMatrix_View cuspMatrix;
  cuspMatrix C;



  cusp::multiply(A,B,C);

  std::cout << " C.column_indices.size():" <<  C.column_indices.size() << std::endl;
  std::cout << " C.values.size():" <<  C.values.size() << std::endl;
  row_mapC = typename in_row_index_view_type::non_const_type("rowmapC", m + 1);
  entriesC = typename in_nonzero_index_view_type::non_const_type ("EntriesC" ,  C.column_indices.size());
  valuesC = typename in_nonzero_value_view_type::non_const_type ("valuesC" ,  C.values.size());

  Kokkos::parallel_for (my_exec_space (0, m + 1) , CopyArrayToCuspArray<typename in_row_index_view_type::non_const_type,
      idx >(row_mapC, (idx *) thrust::raw_pointer_cast(C.row_offsets.data())));
  Kokkos::parallel_for (my_exec_space (0, C.column_indices.size()) , CopyArrayToCuspArray<typename in_nonzero_index_view_type::non_const_type,
      idx >(entriesC, (idx *) thrust::raw_pointer_cast(C.column_indices.data())));
  Kokkos::parallel_for (my_exec_space (0, C.values.size()) , CopyArrayToCuspArray<typename in_nonzero_value_view_type::non_const_type,
      value_type>(valuesC, (value_type *) thrust::raw_pointer_cast(C.values.data())));

#else
  std::cerr << "CUSP IS NOT DEFINED" << std::endl;
  return;
#endif
}
}
}
}

}
#endif
