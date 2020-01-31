/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef _KOKKOSSPGEMMCUSP_HPP
#define _KOKKOSSPGEMMCUSP_HPP

#ifdef KERNELS_HAVE_CUSP
#include <cusp/multiply.h>
#include <cusp/csr_matrix.h>
#endif

namespace KokkosSparse{
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
  typename ain_row_index_view_type,
  typename ain_nonzero_index_view_type,
  typename ain_nonzero_value_view_type,
  typename bin_row_index_view_type,
  typename bin_nonzero_index_view_type,
  typename bin_nonzero_value_view_type,
  typename cin_row_index_view_type,
  typename cin_nonzero_index_view_type,
  typename cin_nonzero_value_view_type>
void CUSP_apply(
    KernelHandle *handle,
    typename KernelHandle::nnz_lno_t m,
    typename KernelHandle::nnz_lno_t n,
    typename KernelHandle::nnz_lno_t k,
    ain_row_index_view_type row_mapA,
    ain_nonzero_index_view_type entriesA,
    ain_nonzero_value_view_type valuesA,

    bool /* transposeA */,
    bin_row_index_view_type row_mapB,
    bin_nonzero_index_view_type entriesB,
    bin_nonzero_value_view_type valuesB,
    bool /* transposeB */,
    cin_row_index_view_type row_mapC,
    cin_nonzero_index_view_type &entriesC,
    cin_nonzero_value_view_type &valuesC){
#ifdef KERNELS_HAVE_CUSP
  typedef typename KernelHandle::nnz_lno_t idx;
  typedef typename KernelHandle::nnz_scalar_t value_type;

  typedef typename ain_row_index_view_type::device_type device1;
  typedef typename ain_nonzero_index_view_type::device_type device2;
  typedef typename ain_nonzero_value_view_type::device_type device3;

  if (Kokkos::Impl::is_same<Kokkos::Cuda, device1 >::value){
    throw std::runtime_error ("MEMORY IS NOT ALLOCATED IN GPU DEVICE for CUSP\n");
    //return;
  }
  if (Kokkos::Impl::is_same<Kokkos::Cuda, device2 >::value){
    throw std::runtime_error ("MEMORY IS NOT ALLOCATED IN GPU DEVICE for CUSP\n");
    //return;
  }
  if (Kokkos::Impl::is_same<Kokkos::Cuda, device3 >::value){
    throw std::runtime_error ("MEMORY IS NOT ALLOCATED IN GPU DEVICE for CUSP\n");
    //return;
  }

  //typedef in_row_index_view_type idx_array_type;

  typedef typename Kokkos::RangePolicy<typename KernelHandle::HandleExecSpace> my_exec_space;




  idx nnzA = entriesA.extent(0);
  idx nnzB = entriesB.extent(0);


  idx *a_xadj = (idx *)row_mapA.data();
  idx *b_xadj = (idx *)row_mapB.data();


  idx *a_adj = (idx *)entriesA.data();
  idx *b_adj = (idx *)entriesB.data();

  value_type *a_ew = valuesA.data();
  value_type *b_ew = valuesB.data();

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

  cuspMatrix_View A(m, n, entriesA.extent(0), arraya_xadj, arraya_adj, arraya_ew);
  cuspMatrix_View B(n, k, entriesB.extent(0), arrayb_xadj, arrayb_adj, arrayb_ew);

  /*
  CopyArrayToCuspArray<typename cuspMatrix::row_offsets_array_type, typename KernelHandle::idx_array_type> Aforward(A.row_offsets, row_mapA);
  Kokkos::parallel_for (my_exec_space (0, m + 1) , Aforward);
  Kokkos::parallel_for (my_exec_space (0, n + 1) , CopyArrayToCuspArray<typename cuspMatrix::row_offsets_array_type, typename KernelHandle::idx_array_type>(B.row_offsets, row_mapB));

  Kokkos::parallel_for (my_exec_space (0, entriesA.extent(0)) , CopyArrayToCuspArray<typename cuspMatrix::column_indices_array_type, typename KernelHandle::idx_edge_array_type>(A.column_indices, entriesA));
  Kokkos::parallel_for (my_exec_space (0, entriesB.extent(0)) , CopyArrayToCuspArray<typename cuspMatrix::column_indices_array_type, typename KernelHandle::idx_edge_array_type>(B.column_indices, entriesB));

  Kokkos::parallel_for (my_exec_space (0, valuesA.extent(0)) , CopyArrayToCuspArray<typename cuspMatrix::values_array_type, typename KernelHandle::value_array_type>(A.values, valuesA));
  Kokkos::parallel_for (my_exec_space (0, valuesB.extent(0)) , CopyArrayToCuspArray<typename cuspMatrix::values_array_type, typename KernelHandle::value_array_type>(B.values, valuesB));
  */

  typedef typename cusp::csr_matrix<idx,value_type,cusp::device_memory> cuspMatrix;
  //typedef cuspMatrix_View cuspMatrix;
  cuspMatrix C;


  Kokkos::Impl::Timer timer1;
  cusp::multiply(A,B,C);
  KernelHandle::HandleExecSpace().fence();
  std::cout << "Actual CUSP SPMM Time:" << timer1.seconds() << std::endl;




  //std::cout << " C.column_indices.size():" <<  C.column_indices.size() << std::endl;
  //std::cout << " C.values.size():" <<  C.values.size() << std::endl;
  //row_mapC = typename cin_row_index_view_type::non_const_type("rowmapC", m + 1);


  handle->set_c_nnz( C.values.size());

  entriesC = typename cin_nonzero_index_view_type::non_const_type (Kokkos::ViewAllocateWithoutInitializing("EntriesC") ,  C.column_indices.size());
  valuesC = typename cin_nonzero_value_view_type::non_const_type (Kokkos::ViewAllocateWithoutInitializing("valuesC"),  C.values.size());

  Kokkos::parallel_for (my_exec_space (0, m + 1) , CopyArrayToCuspArray<typename cin_row_index_view_type::non_const_type,
      idx >(row_mapC, (idx *) thrust::raw_pointer_cast(C.row_offsets.data())));
  Kokkos::parallel_for (my_exec_space (0, C.column_indices.size()) , CopyArrayToCuspArray<typename cin_nonzero_index_view_type::non_const_type,
      idx >(entriesC, (idx *) thrust::raw_pointer_cast(C.column_indices.data())));
  Kokkos::parallel_for (my_exec_space (0, C.values.size()) , CopyArrayToCuspArray<typename cin_nonzero_value_view_type::non_const_type,
      value_type>(valuesC, (value_type *) thrust::raw_pointer_cast(C.values.data())));

#else
  (void)handle;
  (void)m;        (void)n;        (void)k;
  (void)row_mapA; (void)row_mapB; (void)row_mapC;
  (void)entriesA; (void)entriesB; (void)entriesC;
  (void)valuesA;  (void)valuesB;  (void)valuesC;
  throw std::runtime_error ("CUSP IS NOT DEFINED\n");
  //return;
#endif
}
}

}
#endif
