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
    if (entriesB.dimension_0() > max_integer|| entriesA.dimension_0() > max_integer){
      throw std::runtime_error ("MKL requires integer values for size type for SPGEMM. Copying to integer will cause overflow.\n");
      return;
    }



    int *a_adj = (int *)entriesA.ptr_on_device();
    int *b_adj = (int *)entriesB.ptr_on_device();

    int *a_xadj = (int *)row_mapA.ptr_on_device();
    int *b_xadj = (int *)row_mapB.ptr_on_device();
    int *c_xadj = (int *)row_mapC.ptr_on_device();

    if (handle->mkl_convert_to_1base)
    {
      handle->persistent_a_xadj = int_persistent_work_view_t("tmpa", m + 1);
      handle->persistent_b_xadj = int_persistent_work_view_t("tmpb", n + 1);
      handle->persistent_c_xadj = int_persistent_work_view_t("tmpc", m + 1);
      int_persistent_work_view_t a_plus_one ("a_plus_one", entriesA.dimension_0());
      int_persistent_work_view_t b_plus_one ("b_plus_one", entriesB.dimension_0());
      handle->persistent_a_adj = a_plus_one;
      handle->persistent_b_adj = b_plus_one;

      KokkosKernels::Impl::kk_a_times_x_plus_b< int_persistent_work_view_t, in_row_index_view_type,   int, int, MyExecSpace>(m + 1,  handle->persistent_a_xadj, row_mapA,  1, 1);
      KokkosKernels::Impl::kk_a_times_x_plus_b< int_persistent_work_view_t, bin_row_index_view_type,   int, int, MyExecSpace>(n + 1, handle->persistent_b_xadj, row_mapB,  1, 1);
      KokkosKernels::Impl::kk_a_times_x_plus_b<   int_persistent_work_view_t, in_nonzero_index_view_type, int, int, MyExecSpace>(entriesA.dimension_0(), a_plus_one, entriesA,  1, 1);
      KokkosKernels::Impl::kk_a_times_x_plus_b< int_persistent_work_view_t, bin_nonzero_index_view_type,  int, int, MyExecSpace>(entriesB.dimension_0(), b_plus_one, entriesB,  1, 1);


      a_adj = (int *)handle->persistent_a_adj.ptr_on_device();
      b_adj = (int *)handle->persistent_b_adj.ptr_on_device();
      a_xadj = handle->persistent_a_xadj.ptr_on_device();
      b_xadj = handle->persistent_b_xadj.ptr_on_device();
      c_xadj = handle->persistent_c_xadj.ptr_on_device();
    }

    
    
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
    if (verbose)
      std::cout << "Sort:" << sort << " Actual MKL2 Symbolic Time:" << timer1.seconds() << std::endl;

    if (handle->mkl_convert_to_1base){
      KokkosKernels::Impl::kk_a_times_x_plus_b< cin_row_index_view_type, int_persistent_work_view_t,  int, int, MyExecSpace>(m + 1, row_mapC, handle->persistent_c_xadj,  1, -1);
      handle->set_c_nnz(row_mapC(m));
    }
    else {
      handle->set_c_nnz(row_mapC(m) - 1);
    }
  }
  else {
    throw std::runtime_error ("MKL requires local ordinals to be integer.\n");
    return;
  }
#else
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


      int *a_xadj = (int *)row_mapA.ptr_on_device();
      int *b_xadj = (int *)row_mapB.ptr_on_device();
      int *c_xadj = (int *)row_mapC.ptr_on_device();

      int *a_adj = (int *)entriesA.ptr_on_device();
      int *b_adj = (int *)entriesB.ptr_on_device();
      

      const value_type *a_ew = valuesA.ptr_on_device();
      const value_type *b_ew = valuesB.ptr_on_device();

      if (handle->mkl_convert_to_1base)
      {
        int_persistent_work_view_t a_xadj_v, b_xadj_v, c_xadj_v;
        a_xadj = (int *) handle->persistent_a_xadj.ptr_on_device();
        b_xadj = (int *) handle->persistent_b_xadj.ptr_on_device();
        c_xadj = (int *) handle->persistent_c_xadj.ptr_on_device();
        int_persistent_work_view_t a_plus_one =  handle->persistent_a_adj;
        int_persistent_work_view_t b_plus_one =  handle->persistent_b_adj;

        a_adj = (int *)a_plus_one.ptr_on_device();
        b_adj = (int *)b_plus_one.ptr_on_device();
      }


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
      std::cout << "c:" << "entriesC:" << entriesC.dimension_0() << std::endl;
      KokkosKernels::Impl::print_1Dview(row_mapC);
*/
      Kokkos::Impl::Timer timer1;

      if (Kokkos::Impl::is_same<value_type, float>::value){

        mkl_scsrmultcsr(&trans, &request, &sort, &mklm, &mkln, &mklk,
                      (float *)a_ew, a_adj, a_xadj,
                      (float *)b_ew, b_adj, b_xadj,
                      (float *)valuesC.ptr_on_device(), entriesC.ptr_on_device(), c_xadj,
                      &nzmax, &info
                      );
        mkl_free_buffers();
      }
      else if (Kokkos::Impl::is_same<value_type, double>::value){

        mkl_dcsrmultcsr(&trans, &request, &sort, &mklm, &mkln, &mklk,
                      (double *)a_ew, a_adj, a_xadj,
                      (double *)b_ew, b_adj, b_xadj,
                      (double *)valuesC.ptr_on_device(), entriesC.ptr_on_device(), c_xadj,
                      &nzmax, &info
                      );
        mkl_free_buffers();
      }
      else {
        throw std::runtime_error ("MKL requires float or double values. Complex values are not implemented yet.\n");
        return;
      }
      if (verbose)
              std::cout << "Sort:" << sort << " Actual MKL2 Numeric Time:" << timer1.seconds() << std::endl;


      if (handle->mkl_convert_to_1base)
      {
        KokkosKernels::Impl::kk_a_times_x_plus_b< cin_nonzero_index_view_type, cin_nonzero_index_view_type,  int, int, MyExecSpace>(entriesC.dimension_0(), entriesC, entriesC,  1, -1);
      }
    }
    else {
      throw std::runtime_error ("MKL requires local ordinals to be integer.\n");
      return;
    }
#else
    throw std::runtime_error ("MKL IS NOT DEFINED\n");
    //return;
#endif
  }
}
}

#endif
