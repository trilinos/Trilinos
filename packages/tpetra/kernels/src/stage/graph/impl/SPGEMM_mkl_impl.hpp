
#ifndef _KOKKOSSPGEMMMKL_HPP
#define _KOKKOSSPGEMMMKL_HPP

//#define KERNELS_HAVE_MKL

#ifdef KERNELS_HAVE_MKL
#include "mkl_spblas.h"
#endif

#include "KokkosKernelsUtils.hpp"
namespace KokkosKernels{

namespace Experimental{

namespace Graph{
namespace Impl{



  template <typename KernelHandle>
  void mkl_apply(
      KernelHandle *handle,
      typename KernelHandle::idx m,
      typename KernelHandle::idx n,
      typename KernelHandle::idx k,
      typename KernelHandle::idx_array_type row_mapA,
      typename KernelHandle::idx_edge_array_type entriesA,
      typename KernelHandle::value_array_type valuesA,

      bool transposeA,
      typename KernelHandle::idx_array_type row_mapB,
      typename KernelHandle::idx_edge_array_type entriesB,
      typename KernelHandle::value_array_type valuesB,
      bool transposeB,
      typename KernelHandle::idx_array_type &row_mapC,
      typename KernelHandle::idx_edge_array_type &entriesC,
      typename KernelHandle::value_array_type &valuesC){

#ifdef KERNELS_HAVE_MKL
    typedef typename KernelHandle::idx idx;
    typedef typename KernelHandle::idx_array_type idx_array_type;

    typedef typename KernelHandle::value_type value_type;


    typedef typename KernelHandle::idx_device_type device1;
    typedef typename KernelHandle::idx_edge_device_type device2;
    typedef typename KernelHandle::value_type_device_type device3;

    std::cout << "RUNNING MKL" << std::endl;

    if (!Kokkos::Impl::is_same<Kokkos::Cuda, device1 >::value){
      std::cerr << "MEMORY IS NOT ALLOCATED IN HOST DEVICE for MKL" << std::endl;
      return;
    }
    if (!Kokkos::Impl::is_same<Kokkos::Cuda, device2 >::value){
      std::cerr << "MEMORY IS NOT ALLOCATED IN HOST DEVICE for MKL" << std::endl;
      return;
    }
    if (!Kokkos::Impl::is_same<Kokkos::Cuda, device3 >::value){
      std::cerr << "MEMORY IS NOT ALLOCATED IN HOST DEVICE for MKL" << std::endl;
      return;
    }


    if (Kokkos::Impl::is_same<idx, int>::value){
      int *a_xadj = (int *)row_mapA.ptr_on_device();
      int *b_xadj = (int *)row_mapB.ptr_on_device();
      int *c_xadj = (int *)row_mapC.ptr_on_device();

      int *a_adj = (int *)entriesA.ptr_on_device();
      int *b_adj = (int *)entriesB.ptr_on_device();
      int *c_adj = (int *)entriesC.ptr_on_device();

      int nnzA = entriesA.dimension_0();
      int nnzB = entriesB.dimension_0();

      value_type *a_ew = valuesA.ptr_on_device();
      value_type *b_ew = valuesB.ptr_on_device();
      value_type *c_ew = valuesC.ptr_on_device();

      sparse_matrix_t *A;
      sparse_matrix_t *B;
      sparse_matrix_t *C;

      if (Kokkos::Impl::is_same<value_type, float>::value){



        if (SPARSE_STATUS_SUCCESS != mkl_sparse_s_create_csr (A, SPARSE_INDEX_BASE_ZERO, m, n, a_xadj, a_xadj + 1, a_adj, a_ew)){
          std::cerr << "CANNOT CREATE mkl_sparse_s_create_csr A" << std::endl;
          return;
        }

        if (SPARSE_STATUS_SUCCESS != mkl_sparse_s_create_csr (B, SPARSE_INDEX_BASE_ZERO, n, k, b_xadj, b_xadj + 1, b_adj, b_ew)){
          std::cerr << "CANNOT CREATE mkl_sparse_s_create_csr B" << std::endl;
          return;
        }


        sparse_operation_t operation;
        if (transposeA && transposeB){
          operation = SPARSE_OPERATION_TRANSPOSE;
        }
        else if (!(transposeA || transposeB)){
          operation = SPARSE_OPERATION_NON_TRANSPOSE;
        }
        else {
          std::cerr << "Ask both to transpose or non transpose for MKL SPGEMM" << std::endl;
          return;
        }

        if (SPARSE_STATUS_SUCCESS != mkl_sparse_spmm (operation, A, B, C)){
          std::cerr << "CANNOT multiply mkl_sparse_spmm " << std::endl;
          return;
        }
        else{
          row_mapC = idx_array_type("rowmapC", m + 1);
          entriesC = typename KernelHandle::idx_edge_array_type ("EntriesC" ,  C.column_indices.size());
          valuesC = typename KernelHandle::value_array_type ("valuesC" ,  C.values.size());

          KokkosKernels::Experimental::Graph::copy_vector<int *, idx_array_type> (m, C.pointerB, row_mapC);
          idx nnz = row_mapC(m) =  C.pointerE[m - 1];

          KokkosKernels::Experimental::Graph::copy_vector<int *, typename KernelHandle::idx_edge_array_type> (nnz, C.columns, entriesC);
          KokkosKernels::Experimental::Graph::copy_vector<float *, typename KernelHandle::value_array_type> (m, C.values, valuesC);
        }


        if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (A)){
          std::cerr << "CANNOT DESTROY mkl_sparse_destroy A" << std::endl;
          return;
        }

        if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (B)){
          std::cerr << "CANNOT DESTROY mkl_sparse_destroy B" << std::endl;
          return;
        }
      }
      else if (Kokkos::Impl::is_same<value_type, double>::value){

        if (SPARSE_STATUS_SUCCESS != mkl_sparse_d_create_csr (A, SPARSE_INDEX_BASE_ZERO, m, n, a_xadj, a_xadj + 1, a_adj, a_ew)){
          std::cerr << "CANNOT CREATE mkl_sparse_d_create_csr A" << std::endl;
          return;
        }

        if (SPARSE_STATUS_SUCCESS != mkl_sparse_d_create_csr (B, SPARSE_INDEX_BASE_ZERO, n, k, b_xadj, b_xadj + 1, b_adj, b_ew)){
          std::cerr << "CANNOT CREATE mkl_sparse_d_create_csr B" << std::endl;
          return;
        }

        sparse_operation_t operation;
        if (transposeA && transposeB){
          operation = SPARSE_OPERATION_TRANSPOSE;
        }
        else if (!(transposeA || transposeB)){
          operation = SPARSE_OPERATION_NON_TRANSPOSE;
        }
        else {
          std::cerr << "Ask both to transpose or non transpose for MKL SPGEMM" << std::endl;
          return;
        }

        if (SPARSE_STATUS_SUCCESS != mkl_sparse_spmm (operation, A, B, C)){
          std::cerr << "CANNOT multiply mkl_sparse_spmm " << std::endl;
          return;
        }
        else{
          row_mapC = idx_array_type("rowmapC", m + 1);
          entriesC = typename KernelHandle::idx_edge_array_type ("EntriesC" ,  C.column_indices.size());
          valuesC = typename KernelHandle::value_array_type ("valuesC" ,  C.values.size());

          KokkosKernels::Experimental::Graph::copy_vector<int *, idx_array_type> (m, C.pointerB, row_mapC);
          idx nnz = row_mapC(m) =  C.pointerE[m - 1];

          KokkosKernels::Experimental::Graph::copy_vector<int *, typename KernelHandle::idx_edge_array_type> (nnz, C.columns, entriesC);
          KokkosKernels::Experimental::Graph::copy_vector<double *, typename KernelHandle::value_array_type> (m, C.values, valuesC);
        }


        if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (A)){
          std::cerr << "CANNOT DESTROY mkl_sparse_destroy A" << std::endl;
          return;
        }

        if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (B)){
          std::cerr << "CANNOT DESTROY mkl_sparse_destroy B" << std::endl;
          return;
        }

      }
      else {
        std::cerr << "CUSPARSE requires float or double values. cuComplex and cuDoubleComplex are not implemented yet." << std::endl;
        return;
      }
    }
    else {
      std::cerr << "MKL requires integer values" << std::endl;
      return;
    }
#else
    std::cerr << "MKL IS NOT DEFINED" << std::endl;
    return;
#endif
  }
}
}
}
}

#endif
