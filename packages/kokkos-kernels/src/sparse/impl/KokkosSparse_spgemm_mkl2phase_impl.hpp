/*
//@HEADER
// ************************************************************************
//
//          KokkosKernels: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef _KOKKOSSPGEMMMKL2_HPP
#define _KOKKOSSPGEMMMKL2_HPP

//#define KOKKOSKERNELS_ENABLE_TPL_MKL


#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#include "mkl.h"
#endif

#include "KokkosKernels_Utils.hpp"
#include <Kokkos_Concepts.hpp>
#include <vector>

namespace KokkosSparse{
namespace Impl{


template <typename KernelHandle,
typename in_row_index_view_type,
typename in_nonzero_index_view_type,
typename bin_row_index_view_type,
typename bin_nonzero_index_view_type,
typename cin_row_index_view_type>
void mkl2phase_symbolic(
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
  
  typedef typename KernelHandle::HandlePersistentMemorySpace HandlePersistentMemorySpace;

  typedef typename Kokkos::View<int *, HandlePersistentMemorySpace> int_persistent_work_view_t;

  typedef typename KernelHandle::HandleExecSpace MyExecSpace;

  if (Kokkos::Impl::is_same<idx, int>::value){

    int_persistent_work_view_t a_xadj_v, b_xadj_v;

    const int max_integer = 2147483647;
    if (entriesB.extent(0) > max_integer|| entriesA.extent(0) > max_integer){
      throw std::runtime_error ("MKL requires integer values for size type for SPGEMM. Copying to integer will cause overflow.\n");
    }


    int *a_adj = (int *)entriesA.data();
    int *b_adj = (int *)entriesB.data();

    int *a_xadj = (int *)row_mapA.data();
    int *b_xadj = (int *)row_mapB.data();
    int *c_xadj = (int *)row_mapC.data();

    if (handle->mkl_convert_to_1base)
    {
      handle->persistent_a_xadj = int_persistent_work_view_t("tmpa", m + 1);
      handle->persistent_b_xadj = int_persistent_work_view_t("tmpb", n + 1);
      handle->persistent_c_xadj = int_persistent_work_view_t("tmpc", m + 1);
      int_persistent_work_view_t a_plus_one ("a_plus_one", entriesA.extent(0));
      int_persistent_work_view_t b_plus_one ("b_plus_one", entriesB.extent(0));
      handle->persistent_a_adj = a_plus_one;
      handle->persistent_b_adj = b_plus_one;

      KokkosKernels::Impl::kk_a_times_x_plus_b< int_persistent_work_view_t, in_row_index_view_type,   int, int, MyExecSpace>(m + 1,  handle->persistent_a_xadj, row_mapA,  1, 1);
      KokkosKernels::Impl::kk_a_times_x_plus_b< int_persistent_work_view_t, bin_row_index_view_type,   int, int, MyExecSpace>(n + 1, handle->persistent_b_xadj, row_mapB,  1, 1);
      KokkosKernels::Impl::kk_a_times_x_plus_b<   int_persistent_work_view_t, in_nonzero_index_view_type, int, int, MyExecSpace>(entriesA.extent(0), a_plus_one, entriesA,  1, 1);
      KokkosKernels::Impl::kk_a_times_x_plus_b< int_persistent_work_view_t, bin_nonzero_index_view_type,  int, int, MyExecSpace>(entriesB.extent(0), b_plus_one, entriesB,  1, 1);


      a_adj = (int *)handle->persistent_a_adj.data();
      b_adj = (int *)handle->persistent_b_adj.data();
      a_xadj = handle->persistent_a_xadj.data();
      b_xadj = handle->persistent_b_xadj.data();
      c_xadj = handle->persistent_c_xadj.data();
    }

#if __INTEL_MKL__ < 2018
    char trans = 'N';
    MKL_INT request = 1;
    MKL_INT sort = handle->get_mkl_sort_option();
    MKL_INT mklm = m, mkln = n, mklk = k;
    MKL_INT info = 0;

    double *mynullptr = NULL;
    int *mynulladj = NULL;
    const int nzmax = 0;

    /*
    KokkosKernels::Impl::print_1Dview(handle->persistent_a_xadj);
    KokkosKernels::Impl::print_1Dview(a_plus_one);
    KokkosKernels::Impl::print_1Dview(handle->persistent_b_xadj);
    KokkosKernels::Impl::print_1Dview(b_plus_one);
    */
    Kokkos::Impl::Timer timer1;

    mkl_dcsrmultcsr(&trans, &request, &sort, &mklm, &mkln, &mklk,
        mynullptr, a_adj, a_xadj,
        mynullptr, b_adj, b_xadj,
        mynullptr, mynulladj, c_xadj,
        &nzmax, &info);

    if (verbose){ 
      std::cout << "Sort:" << sort << " Actual MKL2 Symbolic Time:" << timer1.seconds() << std::endl; 
    }

    if (handle->mkl_convert_to_1base){
      KokkosKernels::Impl::kk_a_times_x_plus_b< cin_row_index_view_type, int_persistent_work_view_t,  int, int, MyExecSpace>(m + 1, row_mapC, handle->persistent_c_xadj,  1, -1);
      handle->set_c_nnz(row_mapC(m));
    }
    else {
      handle->set_c_nnz(row_mapC(m) - 1);
    }
#endif

#if __INTEL_MKL__ == 2018 && __INTEL_MKL_UPDATE__ >= 2
    MKL_INT mklm = m, mkln = n;
    double *mynullptr = NULL;

    sparse_matrix_t A;
    sparse_matrix_t B;
    sparse_matrix_t C;

    // Goal: Set c_xadj (which is the rowptr) from C

    if (handle->mkl_convert_to_1base) { // a*, b* already converted to 1base above...
      if (SPARSE_STATUS_SUCCESS != mkl_sparse_d_create_csr (&A, SPARSE_INDEX_BASE_ONE, mklm, mkln, a_xadj, a_xadj + 1, a_adj, mynullptr)){
        throw std::runtime_error ("CANNOT CREATE mkl_sparse_s_create_csr A matrix\n");
      }
  
      if (SPARSE_STATUS_SUCCESS != mkl_sparse_d_create_csr (&B, SPARSE_INDEX_BASE_ONE, n, k, b_xadj, b_xadj + 1, b_adj, mynullptr)){
        throw std::runtime_error ("CANNOT CREATE mkl_sparse_s_create_csr B matrix\n");
      }
    } else {
      if (SPARSE_STATUS_SUCCESS != mkl_sparse_d_create_csr (&A, SPARSE_INDEX_BASE_ZERO, mklm, mkln, a_xadj, a_xadj + 1, a_adj, mynullptr)){
        throw std::runtime_error ("CANNOT CREATE mkl_sparse_s_create_csr A matrix\n");
      }
  
      if (SPARSE_STATUS_SUCCESS != mkl_sparse_d_create_csr (&B, SPARSE_INDEX_BASE_ZERO, n, k, b_xadj, b_xadj + 1, b_adj, mynullptr)){
        throw std::runtime_error ("CANNOT CREATE mkl_sparse_s_create_csr B matrix\n");
      }
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
    }

    matrix_descr common_mtx_props;
    common_mtx_props.type = SPARSE_MATRIX_TYPE_GENERAL;
    common_mtx_props.mode = SPARSE_FILL_MODE_FULL;
    common_mtx_props.diag = SPARSE_DIAG_NON_UNIT;

    Kokkos::Impl::Timer timer1;
    // options: SPARSE_STAGE_FULL_MULT vs SPARSE_STAGE_NNZ_COUNT then SPARSE_STAGE_FINALIZE_MULT
    bool success = SPARSE_STATUS_SUCCESS != mkl_sparse_sp2m (operation, common_mtx_props, A, operation, common_mtx_props, B, SPARSE_STAGE_NNZ_COUNT, &C); // success is "true" if mkl_sparse_spmm does not return success

    if (verbose){ 
      std::cout << "Actual DOUBLE MKL SPMM Time:" << timer1.seconds() << std::endl; 
    }

    if (success) {
      throw std::runtime_error ("ERROR at SPGEMM multiplication in mkl_sparse_spmm\n");
    } 
    else {

      // Copy sparse_matrix_t C results back to input data structure
      sparse_index_base_t c_indexing;
      MKL_INT c_rows, c_cols, *rows_end, *columns; // use c_xadj as rows_start
      double *values; // should return null

      if (SPARSE_STATUS_SUCCESS !=
          //mkl_sparse_s_export_csr (C, &c_indexing, &c_rows, &c_cols, &rows_start, &rows_end, &columns, &values)) 
          mkl_sparse_d_export_csr (C, &c_indexing, &c_rows, &c_cols, &c_xadj, &rows_end, &columns, &values)) 
      {
        throw std::runtime_error ("ERROR at exporting result matrix in mkl_sparse_spmm\n");
      }

//      if (SPARSE_INDEX_BASE_ZERO != c_indexing){
//        throw std::runtime_error ("C is not zero based indexed\n");
//      }
      if (handle->mkl_convert_to_1base && (c_indexing == SPARSE_INDEX_BASE_ONE)) { // Need to convert back to base0
        KokkosKernels::Impl::kk_a_times_x_plus_b< cin_row_index_view_type, int_persistent_work_view_t,  int, int, MyExecSpace>(m + 1, row_mapC, handle->persistent_c_xadj,  1, -1);
        handle->set_c_nnz(row_mapC(m));
      }
      else {
        handle->set_c_nnz(row_mapC(m) - 1);
      }

    } // end else !success

    // Cleanup...
    if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (A)){
      throw std::runtime_error ("Error at mkl_sparse_destroy A\n");
    }
    if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (B)){
      throw std::runtime_error ("Error at mkl_sparse_destroy B\n");
    }
    if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (C)){
      throw std::runtime_error ("Error at mkl_sparse_destroy C\n");
    }
#elif __INTEL_MKL__ == 2018 && __INTEL_MKL_UPDATE__ < 2
    throw std::runtime_error ("Intel MKL version 18 must have update 2 - use intel/18.2.xyz\n");
#else
    throw std::runtime_error ("Intel MKL versions > 18 are not yet tested/supported\n");
#endif

  }
  else {
    throw std::runtime_error ("MKL requires local ordinals to be integer.\n");
  }
#else
  throw std::runtime_error ("MKL IS NOT DEFINED\n");
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
  void mkl2phase_apply(
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
      cin_nonzero_index_view_type &entriesC,
      cin_nonzero_value_view_type &valuesC,
      bool verbose = false){

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL

    typedef typename KernelHandle::nnz_lno_t idx;

    typedef typename KernelHandle::HandlePersistentMemorySpace HandlePersistentMemorySpace;

    typedef typename Kokkos::View<int *, HandlePersistentMemorySpace> int_persistent_work_view_t;

    typedef typename KernelHandle::nnz_scalar_t value_type;

    typedef typename KernelHandle::HandleExecSpace MyExecSpace;

    if (Kokkos::Impl::is_same<idx, int>::value){

      int *a_xadj = (int *)row_mapA.data();
      int *b_xadj = (int *)row_mapB.data();
      int *c_xadj = (int *)row_mapC.data();

      int *a_adj = (int *)entriesA.data();
      int *b_adj = (int *)entriesB.data();
      

      if (handle->mkl_convert_to_1base)
      {
        int_persistent_work_view_t a_xadj_v, b_xadj_v, c_xadj_v;
        a_xadj = (int *) handle->persistent_a_xadj.data();
        b_xadj = (int *) handle->persistent_b_xadj.data();
        c_xadj = (int *) handle->persistent_c_xadj.data();
        int_persistent_work_view_t a_plus_one =  handle->persistent_a_adj;
        int_persistent_work_view_t b_plus_one =  handle->persistent_b_adj;

        a_adj = (int *)a_plus_one.data();
        b_adj = (int *)b_plus_one.data();
      }

#if __INTEL_MKL__ < 2018
      const value_type *a_ew = valuesA.data();
      const value_type *b_ew = valuesB.data();

      char trans = 'N';
      MKL_INT request = 2;
      MKL_INT sort = handle->get_mkl_sort_option();
      MKL_INT mklm = m, mkln = n, mklk = k;
      MKL_INT info = 0, nzmax = 2147483647;
/*
      KokkosKernels::Impl::print_1Dview(handle->persistent_a_xadj);
      KokkosKernels::Impl::print_1Dview(a_plus_one);
      KokkosKernels::Impl::print_1Dview(handle->persistent_b_xadj);
      KokkosKernels::Impl::print_1Dview(b_plus_one);
      KokkosKernels::Impl::print_1Dview(handle->persistent_c_xadj);
      KokkosKernels::Impl::print_1Dview(valuesA);
      KokkosKernels::Impl::print_1Dview(valuesB);


      std::cout << "A" << std::endl;
      KokkosKernels::Impl::print_1Dview(row_mapA);
      KokkosKernels::Impl::print_1Dview(entriesA);
      std::cout << "B" << std::endl;
      KokkosKernels::Impl::print_1Dview(row_mapB);
      KokkosKernels::Impl::print_1Dview(entriesB);
      std::cout << "c:" << "entriesC:" << entriesC.extent(0) << std::endl;
      KokkosKernels::Impl::print_1Dview(row_mapC);
*/
      Kokkos::Impl::Timer timer1;

      if (Kokkos::Impl::is_same<value_type, float>::value){

        mkl_scsrmultcsr(&trans, &request, &sort, &mklm, &mkln, &mklk,
                      (float *)a_ew, a_adj, a_xadj,
                      (float *)b_ew, b_adj, b_xadj,
                      (float *)valuesC.data(), entriesC.data(), c_xadj,
                      &nzmax, &info
                      );
        mkl_free_buffers();
      }
      else if (Kokkos::Impl::is_same<value_type, double>::value){

        mkl_dcsrmultcsr(&trans, &request, &sort, &mklm, &mkln, &mklk,
                      (double *)a_ew, a_adj, a_xadj,
                      (double *)b_ew, b_adj, b_xadj,
                      (double *)valuesC.data(), entriesC.data(), c_xadj,
                      &nzmax, &info
                      );
        mkl_free_buffers();
      }
      else {
        throw std::runtime_error ("MKL requires float or double values. Complex values are not implemented yet.\n");
      }
      if (verbose)
      { std::cout << "Sort:" << sort << " Actual MKL2 Numeric Time:" << timer1.seconds() << std::endl; }


      if (handle->mkl_convert_to_1base)
      {
        KokkosKernels::Impl::kk_a_times_x_plus_b< cin_nonzero_index_view_type, cin_nonzero_index_view_type,  int, int, MyExecSpace>(entriesC.extent(0), entriesC, entriesC,  1, -1);
      }
#endif

#if __INTEL_MKL__ == 2018 && __INTEL_MKL_UPDATE__ >= 2
      value_type *a_ew = const_cast<value_type*>(valuesA.data());
      value_type *b_ew = const_cast<value_type*>(valuesB.data());

      char trans = 'N';
      MKL_INT request = 2;
      MKL_INT sort = handle->get_mkl_sort_option();
      MKL_INT mklm = m, mkln = n, mklk = k;
      MKL_INT info = 0, nzmax = 2147483647;


      sparse_matrix_t A;
      sparse_matrix_t B;
      sparse_matrix_t C;

      if (handle->mkl_convert_to_1base) { // a*, b* already converted to 1base above...
        if (Kokkos::Impl::is_same<value_type, double>::value){
          if (SPARSE_STATUS_SUCCESS != mkl_sparse_d_create_csr (&A, SPARSE_INDEX_BASE_ONE, mklm, mkln, a_xadj, a_xadj + 1, a_adj, reinterpret_cast<double*>(a_ew))){
            throw std::runtime_error ("CANNOT CREATE mkl_sparse_s_create_csr A matrix\n");
          }
        }
        else if (Kokkos::Impl::is_same<value_type, float>::value){
          if (SPARSE_STATUS_SUCCESS != mkl_sparse_s_create_csr (&A, SPARSE_INDEX_BASE_ONE, mklm, mkln, a_xadj, a_xadj + 1, a_adj, reinterpret_cast<float*>(a_ew))){
            throw std::runtime_error ("CANNOT CREATE mkl_sparse_s_create_csr A matrix\n");
          }
        }
    
        if (Kokkos::Impl::is_same<value_type, double>::value){
          if (SPARSE_STATUS_SUCCESS != mkl_sparse_d_create_csr (&B, SPARSE_INDEX_BASE_ONE, n, k, b_xadj, b_xadj + 1, b_adj, reinterpret_cast<double*>(b_ew))){
            throw std::runtime_error ("CANNOT CREATE mkl_sparse_s_create_csr B matrix\n");
          }
        }
        else if (Kokkos::Impl::is_same<value_type, float>::value){
          if (SPARSE_STATUS_SUCCESS != mkl_sparse_s_create_csr (&B, SPARSE_INDEX_BASE_ONE, n, k, b_xadj, b_xadj + 1, b_adj, reinterpret_cast<float*>(b_ew))){
            throw std::runtime_error ("CANNOT CREATE mkl_sparse_s_create_csr B matrix\n");
          }
        }
      } else {
        if (Kokkos::Impl::is_same<value_type, double>::value){
          if (SPARSE_STATUS_SUCCESS != mkl_sparse_d_create_csr (&A, SPARSE_INDEX_BASE_ZERO, mklm, mkln, a_xadj, a_xadj + 1, a_adj, reinterpret_cast<double*>(a_ew))){
            throw std::runtime_error ("CANNOT CREATE mkl_sparse_s_create_csr A matrix\n");
          }
        }
        else if (Kokkos::Impl::is_same<value_type, float>::value){
          if (SPARSE_STATUS_SUCCESS != mkl_sparse_s_create_csr (&A, SPARSE_INDEX_BASE_ZERO, mklm, mkln, a_xadj, a_xadj + 1, a_adj, reinterpret_cast<float*>(a_ew))){
            throw std::runtime_error ("CANNOT CREATE mkl_sparse_s_create_csr A matrix\n");
          }
        }
    
        if (Kokkos::Impl::is_same<value_type, double>::value){
          if (SPARSE_STATUS_SUCCESS != mkl_sparse_d_create_csr (&B, SPARSE_INDEX_BASE_ZERO, n, k, b_xadj, b_xadj + 1, b_adj, reinterpret_cast<double*>(b_ew))){
            throw std::runtime_error ("CANNOT CREATE mkl_sparse_s_create_csr B matrix\n");
          }
        }
        else if (Kokkos::Impl::is_same<value_type, float>::value){
          if (SPARSE_STATUS_SUCCESS != mkl_sparse_s_create_csr (&B, SPARSE_INDEX_BASE_ZERO, n, k, b_xadj, b_xadj + 1, b_adj, reinterpret_cast<float*>(b_ew))){
            throw std::runtime_error ("CANNOT CREATE mkl_sparse_s_create_csr B matrix\n");
          }
        }
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
      }
  
      matrix_descr common_mtx_props;
      common_mtx_props.type = SPARSE_MATRIX_TYPE_GENERAL;
      common_mtx_props.mode = SPARSE_FILL_MODE_FULL;
      common_mtx_props.diag = SPARSE_DIAG_NON_UNIT;

      Kokkos::Impl::Timer timer1;
      // options: SPARSE_STAGE_FULL_MULT vs SPARSE_STAGE_NNZ_COUNT then SPARSE_STAGE_FINALIZE_MULT
      bool success = SPARSE_STATUS_SUCCESS != mkl_sparse_sp2m (operation, common_mtx_props, A, operation, common_mtx_props, B, SPARSE_STAGE_FINALIZE_MULT, &C); // success is "true" if mkl_sparse_spmm does not return success

      if (verbose){ 
        std::cout << "Actual MKL SPMM Time:" << timer1.seconds() << std::endl; 
      }

      if (success) {
        throw std::runtime_error ("ERROR at SPGEMM multiplication in mkl_sparse_spmm\n");
      }

      // Copy C components back
      sparse_index_base_t c_indexing;
      MKL_INT c_rows, c_cols, *rows_start, *rows_end, *columns;
      typedef double values_type;
      values_type *values;

      if (SPARSE_STATUS_SUCCESS !=
          //mkl_sparse_s_export_csr (C,
          mkl_sparse_d_export_csr (C,
              &c_indexing, &c_rows, &c_cols, &rows_start, &rows_end, &columns, &values)){
        throw std::runtime_error ("ERROR at exporting result matrix in mkl_sparse_spmm\n");
      }

      if (SPARSE_INDEX_BASE_ZERO != c_indexing){
        throw std::runtime_error ("C is not zero based indexed\n");
      }

      //KokkosKernels::Impl::copy_vector<MKL_INT *, typename cin_row_index_view_type::non_const_type, MyExecSpace> (m, rows_start, row_mapC);
      idx nnz = row_mapC(m) =  rows_end[m - 1];

      KokkosKernels::Impl::copy_vector<MKL_INT *, typename cin_nonzero_index_view_type::non_const_type , MyExecSpace> (nnz, columns, entriesC);
      KokkosKernels::Impl::copy_vector<values_type *, typename cin_nonzero_value_view_type::non_const_type, MyExecSpace> (nnz, values, valuesC);


      if (handle->mkl_convert_to_1base)
      {
        KokkosKernels::Impl::kk_a_times_x_plus_b< cin_nonzero_index_view_type, cin_nonzero_index_view_type,  int, int, MyExecSpace>(entriesC.extent(0), entriesC, entriesC,  1, -1);
      }


      // Cleanup
      if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (A)){
        throw std::runtime_error ("Error at mkl_sparse_destroy A\n");
      }
      if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (B)){
        throw std::runtime_error ("Error at mkl_sparse_destroy B\n");
      }
      if (SPARSE_STATUS_SUCCESS != mkl_sparse_destroy (C)){
        throw std::runtime_error ("Error at mkl_sparse_destroy C\n");
      }
#elif __INTEL_MKL__ == 2018 && __INTEL_MKL_UPDATE__ < 2
      throw std::runtime_error ("Intel MKL version 18 must have update 2 - use intel/18.2.xyz\n");
#else
      throw std::runtime_error ("Intel MKL versions > 18 are not yet tested/supported\n");
#endif

    }
    else {
      throw std::runtime_error ("MKL requires local ordinals to be integer.\n");
    }
#else
    throw std::runtime_error ("MKL IS NOT DEFINED\n");
#endif
  } // end mkl2phase_apply
} } // namespace KokkosKernels::Impl

#endif
