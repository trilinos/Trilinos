
#ifndef _KOKKOSSPGEMMMKL_HPP
#define _KOKKOSSPGEMMMKL_HPP

//#define KERNELS_HAVE_MKL


#ifdef KERNELS_HAVE_MKL
#include "mkl_spblas.h"
#endif

#include "KokkosKernels_Utils.hpp"
namespace KokkosKernels{

namespace Experimental{

namespace Graph{
namespace Impl{




  template <typename KernelHandle,
  typename in_row_index_view_type,
  typename in_nonzero_index_view_type,
  typename in_nonzero_value_view_type>
  void mkl_apply(
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

#ifdef KERNELS_HAVE_MKL

    typedef typename KernelHandle::row_lno_t idx;
    typedef in_row_index_view_type idx_array_type;

    typedef typename KernelHandle::nnz_scalar_t value_type;


    typedef typename in_row_index_view_type::device_type device1;
    typedef typename in_nonzero_index_view_type::device_type device2;
    typedef typename in_nonzero_value_view_type::device_type device3;

    typedef typename KernelHandle::HandleExecSpace MyExecSpace;

    std::cout << "RUNNING MKL" << std::endl;

#if defined( KOKKOS_HAVE_CUDA )
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
#endif

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

      sparse_matrix_t A;
      sparse_matrix_t B;
      sparse_matrix_t C;

      if (Kokkos::Impl::is_same<value_type, float>::value){



        if (SPARSE_STATUS_SUCCESS != mkl_sparse_s_create_csr (&A, SPARSE_INDEX_BASE_ZERO, m, n, a_xadj, a_xadj + 1, a_adj, (float *)a_ew)){
          std::cerr << "CANNOT CREATE mkl_sparse_s_create_csr A" << std::endl;
          return;
        }

        if (SPARSE_STATUS_SUCCESS != mkl_sparse_s_create_csr (&B, SPARSE_INDEX_BASE_ZERO, n, k, b_xadj, b_xadj + 1, b_adj, (float *)b_ew)){
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


        Kokkos::Impl::Timer timer1;
        bool success = SPARSE_STATUS_SUCCESS != mkl_sparse_spmm (operation, A, B, &C);
        std::cout << "Actual FLOAT MKL SPMM Time:" << timer1.seconds() << std::endl;

        if (success){
          std::cerr << "CANNOT multiply mkl_sparse_spmm " << std::endl;
          return;
        }
        else{

          sparse_index_base_t c_indexing;
          MKL_INT c_rows, c_cols, *rows_start, *rows_end, *columns;
          float *values;

          if (SPARSE_STATUS_SUCCESS !=
              mkl_sparse_s_export_csr (C,
                  &c_indexing, &c_rows, &c_cols, &rows_start, &rows_end, &columns, &values)){
            std::cerr << "CANNOT export result matrix " << std::endl;
            return;
          }

          if (SPARSE_INDEX_BASE_ZERO != c_indexing){
            std::cerr << "C is not zero based indexed." << std::endl;
            return;
          }


          row_mapC = typename in_row_index_view_type::non_const_type(Kokkos::ViewAllocateWithoutInitializing("rowmapC"), c_rows + 1);
          entriesC = typename in_nonzero_index_view_type::non_const_type (Kokkos::ViewAllocateWithoutInitializing("EntriesC") , rows_end[m - 1] );
          valuesC = typename in_nonzero_value_view_type::non_const_type (Kokkos::ViewAllocateWithoutInitializing("valuesC") ,  rows_end[m - 1]);

          KokkosKernels::Experimental::Util::copy_vector<MKL_INT *, typename in_row_index_view_type::non_const_type, MyExecSpace> (m, rows_start, row_mapC);
          idx nnz = row_mapC(m) =  rows_end[m - 1];

          KokkosKernels::Experimental::Util::copy_vector<MKL_INT *, typename in_nonzero_index_view_type::non_const_type , MyExecSpace> (nnz, columns, entriesC);
          KokkosKernels::Experimental::Util::copy_vector<float *, typename in_nonzero_value_view_type::non_const_type, MyExecSpace> (m, values, valuesC);
        }


        if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (A)){
          std::cerr << "CANNOT DESTROY mkl_sparse_destroy A" << std::endl;
          return;
        }

        if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (B)){
          std::cerr << "CANNOT DESTROY mkl_sparse_destroy B" << std::endl;
          return;
        }
        if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (C)){
          std::cerr << "CANNOT DESTROY mkl_sparse_destroy C" << std::endl;
          return;
        }
      }
      else if (Kokkos::Impl::is_same<value_type, double>::value){

        /*
        std::cout << "create a" << std::endl;
        std::cout << "m:" << m << " n:" << n << std::endl;
        std::cout << "a_xadj[0]:" << a_xadj[0] << " a_xadj[m]:" << a_xadj[m] << std::endl;
        std::cout << "a_adj[a_xadj[m] - 1]:" << a_adj[a_xadj[m] - 1] << " a_ew[a_xadj[m] - 1]:" << a_ew[a_xadj[m] - 1] << std::endl;
        */
        if (SPARSE_STATUS_SUCCESS != mkl_sparse_d_create_csr (&A, SPARSE_INDEX_BASE_ZERO, m, n, a_xadj, a_xadj + 1, a_adj, (double *)a_ew)){
          std::cerr << "CANNOT CREATE mkl_sparse_d_create_csr A" << std::endl;
          return;
        }

        //std::cout << "create b" << std::endl;
        if (SPARSE_STATUS_SUCCESS != mkl_sparse_d_create_csr (&B, SPARSE_INDEX_BASE_ZERO, n, k, b_xadj, b_xadj + 1, b_adj, (double *) b_ew)){
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


        Kokkos::Impl::Timer timer1;
        bool success = SPARSE_STATUS_SUCCESS != mkl_sparse_spmm (operation, A, B, &C);
        std::cout << "Actual DOUBLE MKL SPMM Time:" << timer1.seconds() << std::endl;

        if (success){
          std::cerr << "CANNOT multiply mkl_sparse_spmm " << std::endl;
          return;
        }
        else{


          sparse_index_base_t c_indexing;
          MKL_INT c_rows, c_cols, *rows_start, *rows_end, *columns;
          double *values;

          if (SPARSE_STATUS_SUCCESS !=
              mkl_sparse_d_export_csr (C,
                  &c_indexing, &c_rows, &c_cols, &rows_start, &rows_end, &columns, &values)){
            std::cerr << "CANNOT export result matrix " << std::endl;
            return;
          }

          if (SPARSE_INDEX_BASE_ZERO != c_indexing){
            std::cerr << "C is not zero based indexed." << std::endl;
            return;
          }
          {
            Kokkos::Impl::Timer copy_time;
            row_mapC = typename in_row_index_view_type::non_const_type(Kokkos::ViewAllocateWithoutInitializing("rowmapC"), c_rows + 1);
            entriesC = typename in_nonzero_index_view_type::non_const_type (Kokkos::ViewAllocateWithoutInitializing("EntriesC") , rows_end[m - 1] );
            valuesC = typename in_nonzero_value_view_type::non_const_type (Kokkos::ViewAllocateWithoutInitializing("valuesC") ,  rows_end[m - 1]);

            KokkosKernels::Experimental::Util::copy_vector<MKL_INT *, typename in_row_index_view_type::non_const_type, MyExecSpace> (m, rows_start, row_mapC);
            idx nnz = row_mapC(m) =  rows_end[m - 1];

            KokkosKernels::Experimental::Util::copy_vector<MKL_INT *, typename in_nonzero_index_view_type::non_const_type, MyExecSpace> (nnz, columns, entriesC);
            KokkosKernels::Experimental::Util::copy_vector<double *, typename in_nonzero_value_view_type::non_const_type, MyExecSpace> (m, values, valuesC);
            double copy_time_d = copy_time.seconds();
            std::cout << "MKL COPYTIME:" << copy_time_d << std::endl;
          }

        }


        if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (A)){
          std::cerr << "CANNOT DESTROY mkl_sparse_destroy A" << std::endl;
          return;
        }

        if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (B)){
          std::cerr << "CANNOT DESTROY mkl_sparse_destroy B" << std::endl;
          return;
        }
        if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (C)){
          std::cerr << "CANNOT DESTROY mkl_sparse_destroy C" << std::endl;
          return;
        }

      }
      else {
        std::cerr << "CUSPARSE requires float or double values. cuComplex and cuDoubleComplex are not implemented yet." << std::endl;
        return;
      }
    }
    else {

      //int *a_xadj = row_mapA.ptr_on_device();
      std::cerr << "MKL requires integer values" << std::endl;

      if (Kokkos::Impl::is_same<idx, unsigned int>::value){
        std::cerr << "MKL is given unsigned integer" << std::endl;
      }
      else if (Kokkos::Impl::is_same<idx, long>::value){
        std::cerr << "MKL is given long" << std::endl;
      }
      else if (Kokkos::Impl::is_same<idx, const int>::value){
        std::cerr << "MKL is given const int" << std::endl;
      }
      else if (Kokkos::Impl::is_same<idx, unsigned long>::value){
        std::cerr << "MKL is given unsigned long" << std::endl;
      }
      else if (Kokkos::Impl::is_same<idx, const unsigned long>::value){
        std::cerr << "MKL is given const unsigned long" << std::endl;
      }
      else{
        std::cerr << "MKL is given something else" << std::endl;
      }
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
