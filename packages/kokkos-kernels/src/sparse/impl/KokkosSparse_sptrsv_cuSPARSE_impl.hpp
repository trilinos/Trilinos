/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
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

#ifndef _KOKKOSSPTRSVCUSPARSE_HPP
#define _KOKKOSSPTRSVCUSPARSE_HPP

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include "cusparse.h"
#endif
namespace KokkosSparse{
namespace Impl{

  template <typename KernelHandle,
  typename ain_row_index_view_type,
  typename ain_nonzero_index_view_type,
  typename ain_values_scalar_view_type>
  void sptrsvcuSPARSE_symbolic(
      KernelHandle *sptrsv_handle,
      typename KernelHandle::nnz_lno_t nrows,
      ain_row_index_view_type row_map,
      ain_nonzero_index_view_type entries,
      ain_values_scalar_view_type values,
      bool trans
      )
  {

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  typedef typename KernelHandle::nnz_lno_t idx_type;
  typedef typename KernelHandle::size_type  size_type;
  typedef typename KernelHandle::scalar_t  scalar_type;
  typedef typename KernelHandle::memory_space  memory_space;

  const bool is_cuda_space = std::is_same<memory_space, Kokkos::CudaSpace>::value || std::is_same<memory_space, Kokkos::CudaUVMSpace>::value || std::is_same<memory_space, Kokkos::CudaHostPinnedSpace>::value;

  if (!is_cuda_space) {
    throw std::runtime_error ("KokkosKernels sptrsvcuSPARSE_symbolic: MEMORY IS NOT ALLOCATED IN GPU DEVICE for CUSPARSE\n");
  }
  else if (std::is_same<idx_type, int>::value) {

    bool is_lower = sptrsv_handle->is_lower_tri();
    sptrsv_handle->create_cuSPARSE_Handle(trans, is_lower);

    typename KernelHandle::SPTRSVcuSparseHandleType *h = sptrsv_handle->get_cuSparseHandle();

    cusparseStatus_t status;
    status = cusparseCreateCsrsv2Info(&(h->info));
    if (CUSPARSE_STATUS_SUCCESS != status)
      std::cout << "csrsv2info create status error name " << (status) << std::endl;

    // query how much memory used in csrsv2, and allocate the buffer
    int nnz = entries.extent_int(0);
    int pBufferSize;

    if (!std::is_same<size_type, int>::value)
      sptrsv_handle->allocate_tmp_int_rowmap(row_map.extent(0));
    const int* rm  = !std::is_same<size_type, int>::value ? sptrsv_handle->get_int_rowmap_ptr_copy(row_map) : (const int*)row_map.data();
    const int* ent =  entries.data();
    const scalar_type* vals = values.data();

    if (std::is_same<scalar_type,double>::value) {
    cusparseDcsrsv2_bufferSize(
      h->handle,
      h->transpose,
      nrows,
      nnz,
      h->descr,
      (double*)vals,
      (int*)rm,
      (int*)ent,
      h->info,
      &pBufferSize);


      // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
      cudaError_t my_error;
      my_error = cudaMalloc((void**)&(h->pBuffer), pBufferSize);

      if (cudaSuccess != my_error)
        std::cout << "cudmalloc pBuffer error_t error name " << cudaGetErrorString(my_error) << std::endl;

      status = cusparseDcsrsv2_analysis(
        h->handle,
        h->transpose,
        nrows,
        nnz,
        h->descr,
        (double*)vals,
        (int*)rm,
        (int*)ent,
        h->info,
        h->policy,
        h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status)
        std::cout << "analysis status error name " << (status) << std::endl;
    }
    else if (std::is_same<scalar_type,float>::value) {
    cusparseScsrsv2_bufferSize(
      h->handle,
      h->transpose,
      nrows,
      nnz,
      h->descr,
      (float*)vals,
      (int*)rm,
      (int*)ent,
      h->info,
      &pBufferSize);


      // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
      cudaError_t my_error;
      my_error = cudaMalloc((void**)&(h->pBuffer), pBufferSize);

      if (cudaSuccess != my_error)
        std::cout << "cudmalloc pBuffer error_t error name " << cudaGetErrorString(my_error) << std::endl;

      status = cusparseScsrsv2_analysis(
        h->handle,
        h->transpose,
        nrows,
        nnz,
        h->descr,
        (float*)vals,
        (int*)rm,
        (int*)ent,
        h->info,
        h->policy,
        h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status)
        std::cout << "analysis status error name " << (status) << std::endl;
    }
    else if (std::is_same<scalar_type,Kokkos::complex<double>>::value) {
    cusparseZcsrsv2_bufferSize(
      h->handle,
      h->transpose,
      nrows,
      nnz,
      h->descr,
      (cuDoubleComplex*)vals,
      (int*)rm,
      (int*)ent,
      h->info,
      &pBufferSize);

      // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
      cudaError_t my_error;
      my_error = cudaMalloc((void**)&(h->pBuffer), pBufferSize);

      if (cudaSuccess != my_error)
        std::cout << "cudmalloc pBuffer error_t error name " << cudaGetErrorString(my_error) << std::endl;

      status = cusparseZcsrsv2_analysis(
        h->handle,
        h->transpose,
        nrows,
        nnz,
        h->descr,
        (cuDoubleComplex*)vals,
        (int*)rm,
        (int*)ent,
        h->info,
        h->policy,
        h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status)
        std::cout << "analysis status error name " << (status) << std::endl;
    }
    else if (std::is_same<scalar_type,Kokkos::complex<float>>::value) {
    cusparseCcsrsv2_bufferSize(
      h->handle,
      h->transpose,
      nrows,
      nnz,
      h->descr,
      (cuComplex*)vals,
      (int*)rm,
      (int*)ent,
      h->info,
      &pBufferSize);

      // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
      cudaError_t my_error;
      my_error = cudaMalloc((void**)&(h->pBuffer), pBufferSize);

      if (cudaSuccess != my_error)
        std::cout << "cudmalloc pBuffer error_t error name " << cudaGetErrorString(my_error) << std::endl;

      status = cusparseCcsrsv2_analysis(
        h->handle,
        h->transpose,
        nrows,
        nnz,
        h->descr,
        (cuComplex*)vals,
        (int*)rm,
        (int*)ent,
        h->info,
        h->policy,
        h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status)
        std::cout << "analysis status error name " << (status) << std::endl;
    }
    else {
      throw std::runtime_error ("CUSPARSE wrapper error: unsupported type.\n");
    }
  }
  else {
    throw std::runtime_error ("CUSPARSE requires local ordinals to be integer.\n");
  }
#else
    (void)sptrsv_handle;
    (void)nrows;
    (void)row_map;
    (void)entries;
    (void)values;
    (void)trans;
    throw std::runtime_error ("CUSPARSE IS NOT DEFINED\n");
    //return;
#endif

  }


  template <typename KernelHandle,
  typename ain_row_index_view_type,
  typename ain_nonzero_index_view_type,
  typename ain_values_scalar_view_type,
  typename b_values_scalar_view_type,
  typename x_values_scalar_view_type>
  void sptrsvcuSPARSE_solve(
      KernelHandle *sptrsv_handle,
      typename KernelHandle::nnz_lno_t nrows,
      ain_row_index_view_type row_map,
      ain_nonzero_index_view_type entries,
      ain_values_scalar_view_type values,
      b_values_scalar_view_type rhs,
      x_values_scalar_view_type lhs,
      bool trans
      )
  {

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  typedef typename KernelHandle::nnz_lno_t idx_type;
  typedef typename KernelHandle::size_type  size_type;
  typedef typename KernelHandle::scalar_t  scalar_type;

  if (std::is_same<idx_type, int>::value) {

    cusparseStatus_t status;

    typename KernelHandle::SPTRSVcuSparseHandleType *h = sptrsv_handle->get_cuSparseHandle();

    int nnz = entries.extent_int(0);

    const int* rm  = !std::is_same<size_type, int>::value ? sptrsv_handle->get_int_rowmap_ptr() : (const int*)row_map.data();
    const int* ent =  entries.data(); 
    const scalar_type* vals = values.data();
    const scalar_type* bv = rhs.data();
    scalar_type* xv = lhs.data();


    if (std::is_same<scalar_type,double>::value) {

      if (h->pBuffer == nullptr) { std::cout << "  pBuffer invalid" << std::endl; }
      const double alpha = double(1);

      status = cusparseDcsrsv2_solve(
        h->handle,
        h->transpose,
        nrows,
        nnz,
        &alpha,
        h->descr,
        (double*)vals,
        (int*)rm,
        (int*)ent,
        h->info,
        (double*)bv,
        (double*)xv,
        h->policy,
        h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status)
        std::cout << "solve status error name " << (status) << std::endl;
    }
    else if (std::is_same<scalar_type,float>::value) {

      if (h->pBuffer == nullptr) { std::cout << "  pBuffer invalid" << std::endl; }
      const float alpha = float(1);

      status = cusparseScsrsv2_solve(
        h->handle,
        h->transpose,
        nrows,
        nnz,
        &alpha,
        h->descr,
        (float*)vals,
        (int*)rm,
        (int*)ent,
        h->info,
        (float*)bv,
        (float*)xv,
        h->policy,
        h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status)
        std::cout << "solve status error name " << (status) << std::endl;
    }
    else if (std::is_same<scalar_type,Kokkos::complex<double>>::value) {
      cuDoubleComplex cualpha;
      cualpha.x = 1.0;
      cualpha.y = 0.0;
      status = cusparseZcsrsv2_solve(
        h->handle,
        h->transpose,
        nrows,
        nnz,
        &cualpha,
        h->descr,
        (cuDoubleComplex*)vals,
        (int*)rm,
        (int*)ent,
        h->info,
        (cuDoubleComplex*)bv,
        (cuDoubleComplex*)xv,
        h->policy,
        h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status)
        std::cout << "solve status error name " << (status) << std::endl;
    }
    else if (std::is_same<scalar_type,Kokkos::complex<float>>::value) {
      cuComplex cualpha;
      cualpha.x = 1.0;
      cualpha.y = 0.0;
      status = cusparseCcsrsv2_solve(
        h->handle,
        h->transpose,
        nrows,
        nnz,
        &cualpha,
        h->descr,
        (cuComplex*)vals,
        (int*)rm,
        (int*)ent,
        h->info,
        (cuComplex*)bv,
        (cuComplex*)xv,
        h->policy,
        h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status)
        std::cout << "solve status error name " << (status) << std::endl;
    }
    else {
      throw std::runtime_error ("CUSPARSE wrapper error: unsupported type.\n");
    }

  }
  else {
    throw std::runtime_error ("CUSPARSE requires local ordinals to be integer.\n");
  }
#else
    (void)sptrsv_handle;
    (void)nrows;
    (void)row_map;
    (void)entries;
    (void)values;
    (void)rhs;
    (void)lhs;
    (void)trans;
    throw std::runtime_error ("CUSPARSE IS NOT DEFINED\n");
#endif

  }

}
}

#endif
