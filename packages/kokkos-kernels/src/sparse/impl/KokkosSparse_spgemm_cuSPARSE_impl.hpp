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

#ifndef _KOKKOSSPGEMMCUSPARSE_HPP
#define _KOKKOSSPGEMMCUSPARSE_HPP

//#define KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include "cusparse.h"
#endif
namespace KokkosSparse{

namespace Impl{


  template <typename KernelHandle,
  typename ain_row_index_view_type,
  typename ain_nonzero_index_view_type,
  typename bin_row_index_view_type,
  typename bin_nonzero_index_view_type,
  typename cin_row_index_view_type>
  void cuSPARSE_symbolic(
      KernelHandle *handle,
      typename KernelHandle::nnz_lno_t m,
      typename KernelHandle::nnz_lno_t n,
      typename KernelHandle::nnz_lno_t k,
      ain_row_index_view_type row_mapA,
      ain_nonzero_index_view_type entriesA,

      bool transposeA,
      bin_row_index_view_type row_mapB,
      bin_nonzero_index_view_type entriesB,
      bool transposeB,
      cin_row_index_view_type row_mapC
      ){

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

    typedef typename ain_row_index_view_type::device_type device1;
    typedef typename ain_nonzero_index_view_type::device_type device2;

    typedef typename KernelHandle::nnz_lno_t idx;
    typedef typename ain_row_index_view_type::non_const_type idx_array_type;


    //TODO this is not correct, check memory space.
    if (Kokkos::Impl::is_same<Kokkos::Cuda, device1 >::value){
      throw std::runtime_error ("MEMORY IS NOT ALLOCATED IN GPU DEVICE for CUSPARSE\n");
      //return;
    }
    if (Kokkos::Impl::is_same<Kokkos::Cuda, device2 >::value){
      throw std::runtime_error ("MEMORY IS NOT ALLOCATED IN GPU DEVICE for CUSPARSE\n");
      //return;
    }

    if (Kokkos::Impl::is_same<idx, int>::value){

      const idx *a_xadj = (int *)row_mapA.ptr_on_device();
      const idx *b_xadj = (int *)row_mapB.ptr_on_device();
      idx *c_xadj = (int *)row_mapC.ptr_on_device();

      const idx *a_adj = entriesA.ptr_on_device();
      const idx *b_adj = entriesB.ptr_on_device();
      handle->create_cuSPARSE_Handle(transposeA, transposeB);
      typename KernelHandle::SPGEMMcuSparseHandleType *h = handle->get_cuSparseHandle();

      int nnzA = entriesA.dimension_0();
      int nnzB = entriesB.dimension_0();

      int baseC, nnzC;
      int *nnzTotalDevHostPtr = &nnzC;

      cusparseXcsrgemmNnz(h->handle,
                          h->transA,
                          h->transB,
                          (int)m,
                          (int)n,
                          (int)k,
                          h->a_descr,
                          nnzA,
                          (int *) a_xadj,
                          (int *)a_adj,
                          h->b_descr,
                          nnzB,
                          (int *)b_xadj,
                          (int *)b_adj,
                          h->c_descr,
                          (int *)c_xadj,
                          nnzTotalDevHostPtr );

      if (NULL != nnzTotalDevHostPtr){
          nnzC = *nnzTotalDevHostPtr;
      }else{
          cudaMemcpy(&nnzC, c_xadj+m, sizeof(int), cudaMemcpyDeviceToHost);
          cudaMemcpy(&baseC, c_xadj, sizeof(int), cudaMemcpyDeviceToHost);
          nnzC -= baseC;
      }
      handle->set_c_nnz(nnzC);
      //entriesC = cin_nonzero_index_view_type(Kokkos::ViewAllocateWithoutInitializing("entriesC"), nnzC);
    }
    else {
      throw std::runtime_error ("CUSPARSE requires local ordinals to be integer.\n");
      //return;
    }
#else
    throw std::runtime_error ("CUSPARSE IS NOT DEFINED\n");
    //return;
#endif

  }



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
  void cuSPARSE_apply(
      KernelHandle *handle,
      typename KernelHandle::nnz_lno_t m,
      typename KernelHandle::nnz_lno_t n,
      typename KernelHandle::nnz_lno_t k,
      ain_row_index_view_type row_mapA,
      ain_nonzero_index_view_type entriesA,
      ain_nonzero_value_view_type valuesA,

      bool transposeA,
      bin_row_index_view_type row_mapB,
      bin_nonzero_index_view_type entriesB,
      bin_nonzero_value_view_type valuesB,
      bool transposeB,
      cin_row_index_view_type row_mapC,
      cin_nonzero_index_view_type entriesC,
      cin_nonzero_value_view_type valuesC){

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
    typedef typename KernelHandle::nnz_lno_t idx;
    typedef ain_row_index_view_type idx_array_type;

    typedef typename KernelHandle::nnz_scalar_t value_type;


    typedef typename ain_row_index_view_type::device_type device1;
    typedef typename ain_nonzero_index_view_type::device_type device2;
    typedef typename ain_nonzero_value_view_type::device_type device3;


    if (Kokkos::Impl::is_same<Kokkos::Cuda, device1 >::value){
      throw std::runtime_error ("MEMORY IS NOT ALLOCATED IN GPU DEVICE for CUSPARSE\n");
      //return;
    }
    if (Kokkos::Impl::is_same<Kokkos::Cuda, device2 >::value){
      throw std::runtime_error ("MEMORY IS NOT ALLOCATED IN GPU DEVICE for CUSPARSE\n");
      //return;
    }
    if (Kokkos::Impl::is_same<Kokkos::Cuda, device3 >::value){
      throw std::runtime_error ("MEMORY IS NOT ALLOCATED IN GPU DEVICE for CUSPARSE\n");
      //return;
    }



    if (Kokkos::Impl::is_same<idx, int>::value){
      int *a_xadj = (int *)row_mapA.ptr_on_device();
      int *b_xadj = (int *)row_mapB.ptr_on_device();
      int *c_xadj = (int *)row_mapC.ptr_on_device();

      int *a_adj = (int *)entriesA.ptr_on_device();
      int *b_adj = (int *)entriesB.ptr_on_device();
      int *c_adj = (int *)entriesC.ptr_on_device();


      typename KernelHandle::SPGEMMcuSparseHandleType *h = handle->get_cuSparseHandle();

      int nnzA = entriesA.dimension_0();
      int nnzB = entriesB.dimension_0();

      value_type *a_ew = (value_type *)valuesA.ptr_on_device();
      value_type *b_ew = (value_type *)valuesB.ptr_on_device();
      value_type *c_ew = (value_type *)valuesC.ptr_on_device();

      if (Kokkos::Impl::is_same<value_type, float>::value){
        cusparseScsrgemm(
            h->handle,
            h->transA,
            h->transB,
            m,
            n,
            k,
            h->a_descr,
            nnzA,
            (float *)a_ew,
            a_xadj,
            a_adj,
            h->b_descr,
            nnzB,
            (float *)b_ew,
            b_xadj,
            b_adj,
            h->c_descr,
            (float *)c_ew,
            c_xadj,
            c_adj);
      }
      else if (Kokkos::Impl::is_same<value_type, double>::value){
        cusparseDcsrgemm(
            h->handle,
            h->transA,
            h->transB,
            m,
            n,
            k,
            h->a_descr,
            nnzA,
            (double *)a_ew,
            a_xadj,
            a_adj,
            h->b_descr,
            nnzB,
            (double *)b_ew,
            b_xadj,
            b_adj,
            h->c_descr,
            (double *)c_ew,
            c_xadj,
            c_adj);
      }
      else {
        throw std::runtime_error ("CUSPARSE requires float or double values. cuComplex and cuDoubleComplex are not implemented yet.\n");
        //return;
      }




    }
    else {
      throw std::runtime_error ("CUSPARSE requires local ordinals to be integer.\n");
      //return;
    }
#else
    throw std::runtime_error ("CUSPARSE IS NOT DEFINED\n");
    //return;
#endif
  }
}
}

#endif
