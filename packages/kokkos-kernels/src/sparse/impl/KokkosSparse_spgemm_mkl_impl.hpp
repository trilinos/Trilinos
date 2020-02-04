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

#ifndef _KOKKOSSPGEMMMKL_HPP
#define _KOKKOSSPGEMMMKL_HPP


#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#include "mkl_spblas.h"
#include "mkl.h"
#endif

#include "KokkosKernels_Utils.hpp"
#include <Kokkos_Concepts.hpp>


namespace KokkosSparse{

namespace Impl{


template <typename KernelHandle,
typename in_row_index_view_type,
typename in_nonzero_index_view_type,
typename bin_row_index_view_type,
typename bin_nonzero_index_view_type,
typename cin_row_index_view_type>
void mkl_symbolic(
    KernelHandle *handle,
    typename KernelHandle::nnz_lno_t m,
    typename KernelHandle::nnz_lno_t n,
    typename KernelHandle::nnz_lno_t k,
    in_row_index_view_type row_mapA,
    in_nonzero_index_view_type entriesA,

    bool transposeA,
    bin_row_index_view_type row_mapB,
    bin_nonzero_index_view_type entriesB,
    bool transposeB,
    cin_row_index_view_type row_mapC,
    bool verbose = false){

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL

  typedef typename KernelHandle::nnz_lno_t idx;
  typedef typename KernelHandle::size_type size_type;


  typedef typename KernelHandle::HandleTempMemorySpace HandleTempMemorySpace;
  typedef typename Kokkos::View<int *, HandleTempMemorySpace> int_temp_work_view_t;



  typedef typename KernelHandle::nnz_scalar_t value_type;






  typedef typename KernelHandle::HandleExecSpace MyExecSpace;
/*
  if (!(
      (Kokkos::Impl::SpaceAccessibility<typename Kokkos::HostSpace::execution_space, typename device1::memory_space>::accessible) &&
      (Kokkos::Impl::SpaceAccessibility<typename Kokkos::HostSpace::execution_space, typename device2::memory_space>::accessible) &&
      (Kokkos::Impl::SpaceAccessibility<typename Kokkos::HostSpace::execution_space, typename device3::memory_space>::accessible) )
      ){
    throw std::runtime_error ("MEMORY IS NOT ALLOCATED IN HOST DEVICE for MKL\n");
    return;
  }
*/
  if (Kokkos::Impl::is_same<idx, int>::value){

    int *a_xadj = NULL;
    int *b_xadj = NULL;
    int_temp_work_view_t a_xadj_v, b_xadj_v;

    if (Kokkos::Impl::is_same<size_type, int>::value){

      a_xadj = (int *)row_mapA.data();
      b_xadj = (int *)row_mapB.data();
    }
    else {


      //TODO test this case.

      Kokkos::Impl::Timer copy_time;
      const int max_integer = 2147483647;
      if (entriesB.extent(0) > max_integer|| entriesA.extent(0) > max_integer){
        throw std::runtime_error ("MKL requires integer values for size type for SPGEMM. Copying to integer will cause overflow.\n");
        return;
      }
      a_xadj_v = int_temp_work_view_t("tmpa", m + 1);
      a_xadj = (int *) a_xadj_v.data();
      b_xadj_v = int_temp_work_view_t("tmpb", n + 1);
      b_xadj = (int *) b_xadj_v.data();

      KokkosKernels::Impl::copy_vector<
          in_row_index_view_type,
          int_temp_work_view_t,
          MyExecSpace> (m+1, row_mapA, a_xadj_v);

      KokkosKernels::Impl::copy_vector<
			bin_row_index_view_type,
          int_temp_work_view_t,
          MyExecSpace> (m+1, row_mapB, b_xadj_v);

      if (verbose)
        std::cout << "MKL COPY size type to int TIME:" << copy_time.seconds() << std::endl;

    }


    int *a_adj = (int *)entriesA.data();
    int *b_adj = (int *)entriesB.data();



    std::vector <value_type> tmp_values (KOKKOSKERNELS_MACRO_MAX(entriesB.extent(0), entriesA.extent(0)));
    value_type *ptmp_values = &(tmp_values[0]);
    value_type *a_ew = ptmp_values;
    value_type *b_ew = ptmp_values;


    sparse_matrix_t A;
    sparse_matrix_t B;
    sparse_matrix_t C;

    if (Kokkos::Impl::is_same<value_type, float>::value){



      if (SPARSE_STATUS_SUCCESS != mkl_sparse_s_create_csr (&A, SPARSE_INDEX_BASE_ZERO, m, n, a_xadj, a_xadj + 1, a_adj, (float *)a_ew)){
        throw std::runtime_error ("CANNOT CREATE mkl_sparse_s_create_csr A matrix\n");
        return;
      }

      if (SPARSE_STATUS_SUCCESS != mkl_sparse_s_create_csr (&B, SPARSE_INDEX_BASE_ZERO, n, k, b_xadj, b_xadj + 1, b_adj, (float *)b_ew)){
        throw std::runtime_error ("CANNOT CREATE mkl_sparse_s_create_csr B matrix\n");
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

        throw std::runtime_error ("MKL either transpose both matrices, or none for SPGEMM\n");
        return;
      }


      Kokkos::Impl::Timer timer1;
      bool success = SPARSE_STATUS_SUCCESS != mkl_sparse_spmm (operation, A, B, &C);
      if (verbose)
      std::cout << "Actual FLOAT MKL SPMM Time in symbolic:" << timer1.seconds() << std::endl;

      if (success){
        throw std::runtime_error ("ERROR at SPGEMM multiplication in mkl_sparse_spmm\n");


        return;
      }
      else{

        sparse_index_base_t c_indexing;
        MKL_INT c_rows, c_cols, *rows_start, *rows_end, *columns;
        float *values;

        if (SPARSE_STATUS_SUCCESS !=
            mkl_sparse_s_export_csr (C,
                &c_indexing, &c_rows, &c_cols, &rows_start, &rows_end, &columns, &values)){
          throw std::runtime_error ("ERROR at exporting result matrix in mkl_sparse_spmm\n");
          return;
        }

        if (SPARSE_INDEX_BASE_ZERO != c_indexing){
          throw std::runtime_error ("C is not zero based indexed\n");
          return;
        }



        KokkosKernels::Impl::copy_vector<MKL_INT *, typename cin_row_index_view_type::non_const_type, MyExecSpace> (m, rows_start, row_mapC);
        idx nnz = row_mapC(m) =  rows_end[m - 1];
        handle->set_c_nnz(nnz);

      }


      if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (A)){
        throw std::runtime_error ("Error at mkl_sparse_destroy A\n");
        return;
      }

      if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (B)){
        throw std::runtime_error ("Error at mkl_sparse_destroy B\n");
        return;
      }
      if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (C)){
        throw std::runtime_error ("Error at mkl_sparse_destroy C\n");
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
        throw std::runtime_error ("CANNOT CREATE mkl_sparse_s_create_csr A matrix\n");
        return;
      }

      //std::cout << "create b" << std::endl;
      if (SPARSE_STATUS_SUCCESS != mkl_sparse_d_create_csr (&B, SPARSE_INDEX_BASE_ZERO, n, k, b_xadj, b_xadj + 1, b_adj, (double *) b_ew)){
        throw std::runtime_error ("CANNOT CREATE mkl_sparse_s_create_csr B matrix\n");
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
        throw std::runtime_error ("MKL either transpose both matrices, or none for SPGEMM\n");
        return;
      }


      Kokkos::Impl::Timer timer1;
      bool success = SPARSE_STATUS_SUCCESS != mkl_sparse_spmm (operation, A, B, &C);
      if (verbose)
      std::cout << "Actual DOUBLE MKL SPMM Time Without Free:" << timer1.seconds() << std::endl;
      mkl_free_buffers();
      if (verbose)
      std::cout << "Actual DOUBLE MKL SPMM Time:" << timer1.seconds() << std::endl;

      if (success){
        throw std::runtime_error ("ERROR at SPGEMM multiplication in mkl_sparse_spmm\n");
        return;
      }
      else{


        sparse_index_base_t c_indexing;
        MKL_INT c_rows, c_cols, *rows_start, *rows_end, *columns;
        double *values;

        if (SPARSE_STATUS_SUCCESS !=
            mkl_sparse_d_export_csr (C,
                &c_indexing, &c_rows, &c_cols, &rows_start, &rows_end, &columns, &values)){
          throw std::runtime_error ("ERROR at exporting result matrix in mkl_sparse_spmm\n");
          return;
        }

        if (SPARSE_INDEX_BASE_ZERO != c_indexing){
          throw std::runtime_error ("C is not zero based indexed\n");
          return;
        }
        if (handle->mkl_keep_output)
        {
          Kokkos::Impl::Timer copy_time;

          KokkosKernels::Impl::copy_vector<MKL_INT *, typename cin_row_index_view_type::non_const_type, MyExecSpace> (m, rows_start, row_mapC);
          idx nnz = row_mapC(m) =  rows_end[m - 1];
          handle->set_c_nnz(nnz);

          double copy_time_d = copy_time.seconds();
          if (verbose)
          std::cout << "MKL COPYTIME:" << copy_time_d << std::endl;
        }

      }


      if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (A)){
        throw std::runtime_error ("Error at mkl_sparse_destroy A\n");
        return;
      }

      if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (B)){
        throw std::runtime_error ("Error at mkl_sparse_destroy B\n");
        return;
      }
      if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (C)){
        throw std::runtime_error ("Error at mkl_sparse_destroy C\n");
        return;
      }

    }
    else {
      throw std::runtime_error ("MKL requires float or double values. Complex values are not implemented yet.\n");
      return;
    }
  }
  else {
    throw std::runtime_error ("MKL requires local ordinals to be integer.\n");
    return;
  }
#else
  (void)handle;
  (void)m;          (void)n;          (void)k;
  (void)row_mapA;   (void)row_mapB;   (void)row_mapC;
  (void)entriesA;   (void)entriesB;
  (void)transposeA; (void)transposeB;
  (void)verbose;
  throw std::runtime_error ("MKL IS NOT DEFINED\n");
  //return;
#endif
}



  template <typename KernelHandle,
  typename in_row_index_view_type,
  typename in_nonzero_index_view_type,
  typename in_nonzero_value_view_type,
  typename bin_row_index_view_type,
  typename bin_nonzero_index_view_type,
  typename bin_nonzero_value_view_type,
  typename cin_row_index_view_type,
  typename cin_nonzero_index_view_type,
  typename cin_nonzero_value_view_type>
  void mkl_apply(
      KernelHandle *handle,
      typename KernelHandle::nnz_lno_t m,
      typename KernelHandle::nnz_lno_t n,
      typename KernelHandle::nnz_lno_t k,
      in_row_index_view_type row_mapA,
      in_nonzero_index_view_type entriesA,
      in_nonzero_value_view_type valuesA,

      bool transposeA,
      bin_row_index_view_type row_mapB,
      bin_nonzero_index_view_type entriesB,
      bin_nonzero_value_view_type valuesB,
      bool transposeB,
      cin_row_index_view_type row_mapC,
      cin_nonzero_index_view_type entriesC,
      cin_nonzero_value_view_type valuesC,
      bool verbose = false){

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL

    typedef typename KernelHandle::nnz_lno_t idx;
    typedef typename KernelHandle::size_type size_type;


    typedef typename KernelHandle::HandleTempMemorySpace HandleTempMemorySpace;
    typedef typename Kokkos::View<int *, HandleTempMemorySpace> int_temp_work_view_t;



    typedef typename KernelHandle::nnz_scalar_t value_type;

    



    typedef typename KernelHandle::HandleExecSpace MyExecSpace;
/*
    if (!(
        (Kokkos::Impl::SpaceAccessibility<typename Kokkos::HostSpace::execution_space, typename device1::memory_space>::accessible) &&
        (Kokkos::Impl::SpaceAccessibility<typename Kokkos::HostSpace::execution_space, typename device2::memory_space>::accessible) &&
        (Kokkos::Impl::SpaceAccessibility<typename Kokkos::HostSpace::execution_space, typename device3::memory_space>::accessible) )
        ){
      throw std::runtime_error ("MEMORY IS NOT ALLOCATED IN HOST DEVICE for MKL\n");
      return;
    }
*/
    if (Kokkos::Impl::is_same<idx, int>::value){

      int *a_xadj = NULL;
      int *b_xadj = NULL;
      int_temp_work_view_t a_xadj_v, b_xadj_v;

      if (Kokkos::Impl::is_same<size_type, int>::value){

        a_xadj = (int *)row_mapA.data();
        b_xadj = (int *)row_mapB.data();
      }
      else {


        //TODO test this case.

        Kokkos::Impl::Timer copy_time;
        const int max_integer = 2147483647;
        if (entriesB.extent(0) > max_integer|| entriesA.extent(0) > max_integer){
          throw std::runtime_error ("MKL requires integer values for size type for SPGEMM. Copying to integer will cause overflow.\n");
          return;
        }
        a_xadj_v = int_temp_work_view_t("tmpa", m + 1);
        a_xadj = (int *) a_xadj_v.data();
        b_xadj_v = int_temp_work_view_t("tmpb", n + 1);
        b_xadj = (int *) b_xadj_v.data();

        KokkosKernels::Impl::copy_vector<
            in_row_index_view_type,
            int_temp_work_view_t,
            MyExecSpace> (m+1, row_mapA, a_xadj_v);

        KokkosKernels::Impl::copy_vector<
			bin_row_index_view_type,
            int_temp_work_view_t,
            MyExecSpace> (m+1, row_mapB, b_xadj_v);

        if (verbose)
          std::cout << "MKL COPY size type to int TIME:" << copy_time.seconds() << std::endl;

      }


      int *a_adj = (int *)entriesA.data();
      int *b_adj = (int *)entriesB.data();


 
      const value_type *a_ew = valuesA.data();
      const value_type *b_ew = valuesB.data();


      sparse_matrix_t A;
      sparse_matrix_t B;
      sparse_matrix_t C;

      if (Kokkos::Impl::is_same<value_type, float>::value){



        if (SPARSE_STATUS_SUCCESS != mkl_sparse_s_create_csr (&A, SPARSE_INDEX_BASE_ZERO, m, n, a_xadj, a_xadj + 1, a_adj, (float *)a_ew)){
          throw std::runtime_error ("CANNOT CREATE mkl_sparse_s_create_csr A matrix\n");
          return;
        }

        if (SPARSE_STATUS_SUCCESS != mkl_sparse_s_create_csr (&B, SPARSE_INDEX_BASE_ZERO, n, k, b_xadj, b_xadj + 1, b_adj, (float *)b_ew)){
          throw std::runtime_error ("CANNOT CREATE mkl_sparse_s_create_csr B matrix\n");
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

          throw std::runtime_error ("MKL either transpose both matrices, or none for SPGEMM\n");
          return;
        }


        Kokkos::Impl::Timer timer1;
        bool success = SPARSE_STATUS_SUCCESS != mkl_sparse_spmm (operation, A, B, &C);
        if (verbose)
        std::cout << "Actual FLOAT MKL SPMM Time:" << timer1.seconds() << std::endl;

        if (success){
          throw std::runtime_error ("ERROR at SPGEMM multiplication in mkl_sparse_spmm\n");


          return;
        }
        else{

          sparse_index_base_t c_indexing;
          MKL_INT c_rows, c_cols, *rows_start, *rows_end, *columns;
          float *values;

          if (SPARSE_STATUS_SUCCESS !=
              mkl_sparse_s_export_csr (C,
                  &c_indexing, &c_rows, &c_cols, &rows_start, &rows_end, &columns, &values)){
            throw std::runtime_error ("ERROR at exporting result matrix in mkl_sparse_spmm\n");
            return;
          }

          if (SPARSE_INDEX_BASE_ZERO != c_indexing){
            throw std::runtime_error ("C is not zero based indexed\n");
            return;
          }


          //KokkosKernels::Impl::copy_vector<MKL_INT *, typename cin_row_index_view_type::non_const_type, MyExecSpace> (m, rows_start, row_mapC);
          idx nnz = row_mapC(m) =  rows_end[m - 1];

          KokkosKernels::Impl::copy_vector<MKL_INT *, typename cin_nonzero_index_view_type::non_const_type , MyExecSpace> (nnz, columns, entriesC);
          KokkosKernels::Impl::copy_vector<float *, typename cin_nonzero_value_view_type::non_const_type, MyExecSpace> (nnz, values, valuesC);
        }


        if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (A)){
          throw std::runtime_error ("Error at mkl_sparse_destroy A\n");
          return;
        }

        if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (B)){
          throw std::runtime_error ("Error at mkl_sparse_destroy B\n");
          return;
        }
        if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (C)){
          throw std::runtime_error ("Error at mkl_sparse_destroy C\n");
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
        if (SPARSE_STATUS_SUCCESS != mkl_sparse_d_create_csr (&A, SPARSE_INDEX_BASE_ZERO, m, n, a_xadj, a_xadj + 1, a_adj, ( double *)a_ew)){
          throw std::runtime_error ("CANNOT CREATE mkl_sparse_s_create_csr A matrix\n");
          return;
        }

        //std::cout << "create b" << std::endl;
        if (SPARSE_STATUS_SUCCESS != mkl_sparse_d_create_csr (&B, SPARSE_INDEX_BASE_ZERO, n, k, b_xadj, b_xadj + 1, b_adj, ( double *) b_ew)){
          throw std::runtime_error ("CANNOT CREATE mkl_sparse_s_create_csr B matrix\n");
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
          throw std::runtime_error ("MKL either transpose both matrices, or none for SPGEMM\n");
          return;
        }


        Kokkos::Impl::Timer timer1;
        bool success = SPARSE_STATUS_SUCCESS != mkl_sparse_spmm (operation, A, B, &C);
        if (verbose)
        std::cout << "Actual DOUBLE MKL SPMM Time Without Free:" << timer1.seconds() << std::endl;

        mkl_free_buffers();
        if (verbose)
        std::cout << "Actual DOUBLE MKL SPMM Time:" << timer1.seconds() << std::endl;

        if (success){
          throw std::runtime_error ("ERROR at SPGEMM multiplication in mkl_sparse_spmm\n");
          return;
        }
        else{


          sparse_index_base_t c_indexing;
          MKL_INT c_rows, c_cols, *rows_start, *rows_end, *columns;
          double *values;

          if (SPARSE_STATUS_SUCCESS !=
              mkl_sparse_d_export_csr (C,
                  &c_indexing, &c_rows, &c_cols, &rows_start, &rows_end, &columns, &values)){
            throw std::runtime_error ("ERROR at exporting result matrix in mkl_sparse_spmm\n");
            return;
          }

          if (SPARSE_INDEX_BASE_ZERO != c_indexing){
            throw std::runtime_error ("C is not zero based indexed\n");
            return;
          }
          if (handle->mkl_keep_output)
          {
            Kokkos::Impl::Timer copy_time;

            //KokkosKernels::Impl::copy_vector<MKL_INT *, typename cin_row_index_view_type::non_const_type, MyExecSpace> (m, rows_start, row_mapC);
            idx nnz = row_mapC(m) =  rows_end[m - 1];

            KokkosKernels::Impl::copy_vector<MKL_INT *, typename cin_nonzero_index_view_type::non_const_type, MyExecSpace> (nnz, columns, entriesC);
            KokkosKernels::Impl::copy_vector<double *, typename cin_nonzero_value_view_type::non_const_type, MyExecSpace> (nnz, values, valuesC);
            double copy_time_d = copy_time.seconds();
            if (verbose)
            std::cout << "MKL COPYTIME:" << copy_time_d << std::endl;
          }

        }


        if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (A)){
          throw std::runtime_error ("Error at mkl_sparse_destroy A\n");
          return;
        }

        if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (B)){
          throw std::runtime_error ("Error at mkl_sparse_destroy B\n");
          return;
        }
        if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (C)){
          throw std::runtime_error ("Error at mkl_sparse_destroy C\n");
          return;
        }

      }
      else {
        throw std::runtime_error ("MKL requires float or double values. Complex values are not implemented yet.\n");
        return;
      }
    }
    else {
      throw std::runtime_error ("MKL requires local ordinals to be integer.\n");
      return;
    }
#else
    (void)handle;
    (void)m;          (void)n;          (void)k;
    (void)row_mapA;   (void)row_mapB;   (void)row_mapC;
    (void)entriesA;   (void)entriesB;   (void)entriesC;
    (void)valuesA;    (void)valuesB;    (void)valuesC;
    (void)transposeA; (void)transposeB;
    (void)verbose;
    throw std::runtime_error ("MKL IS NOT DEFINED\n");
    //return;
#endif
  }
}
}

#endif
