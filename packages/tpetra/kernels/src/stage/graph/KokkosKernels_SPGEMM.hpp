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
#ifndef _KOKKOS_GAUSSSEIDEL_HPP
#define _KOKKOS_GAUSSSEIDEL_HPP

#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_SPGEMM_viennaCL_impl.hpp"
#include "KokkosKernels_SPGEMM_cuSPARSE_impl.hpp"
#include "KokkosKernels_SPGEMM_CUSP_impl.hpp"
#include "KokkosKernels_SPGEMM_mkl_impl.hpp"
#include "KokkosKernels_SPGEMM_impl.hpp"

namespace KokkosKernels{

namespace Experimental{

namespace Graph{

  template <typename KernelHandle,
  typename alno_row_view_t_,
  typename alno_nnz_view_t_,
  typename blno_row_view_t_,
  typename blno_nnz_view_t_,
  typename clno_row_view_t_>
  void spgemm_symbolic(
      KernelHandle *handle,
      typename KernelHandle::nnz_lno_t m,
      typename KernelHandle::nnz_lno_t n,
      typename KernelHandle::nnz_lno_t k,
      alno_row_view_t_ row_mapA,
      alno_nnz_view_t_ entriesA,
      bool transposeA,
      blno_row_view_t_ row_mapB,
      blno_nnz_view_t_ entriesB,
      bool transposeB,
      clno_row_view_t_ row_mapC){

    typedef typename KernelHandle::SPGEMMHandleType spgemmHandleType;
    spgemmHandleType *sh = handle->get_spgemm_handle();
    switch (sh->get_algorithm_type()){

    case SPGEMM_CUSPARSE:
      Impl::cuSPARSE_symbolic
      <spgemmHandleType,
      alno_row_view_t_,
      alno_nnz_view_t_,
      blno_row_view_t_,
      blno_nnz_view_t_,
      clno_row_view_t_>(sh, m,n,k,
          row_mapA, entriesA, transposeA,
          row_mapB, entriesB, transposeB,
          row_mapC);
      break;

    case SPGEMM_CUSP:
      break;

    case SPGEMM_DEFAULT:
    case SPGEMM_KK_MEMSPEED:
    case SPGEMM_KK_SPEED:
    case SPGEMM_KK_MEMORY:
    case SPGEMM_KK_COLOR:
    case SPGEMM_KK_MULTICOLOR:
    case SPGEMM_KK_MULTICOLOR2:
    {
      KokkosKernels::Experimental::Graph::Impl::KokkosSPGEMM
      <KernelHandle,
        alno_row_view_t_, alno_nnz_view_t_, typename KernelHandle::in_scalar_nnz_view_t,
        blno_row_view_t_, blno_nnz_view_t_, typename KernelHandle::in_scalar_nnz_view_t>
      kspgemm (handle,m,n,k,row_mapA, entriesA, transposeA, row_mapB, entriesB, transposeB);
      kspgemm.KokkosSPGEMM_symbolic(row_mapC);
    }
      break;

    case SPGEMM_SERIAL:
    case SPGEMM_DEBUG:
    case SPGEMM_MKL:
    default:
      break;
    }
    sh->set_call_symbolic();

  }


  template <typename KernelHandle,
    typename alno_row_view_t_,
    typename alno_nnz_view_t_,
    typename ascalar_nnz_view_t_,
    typename blno_row_view_t_,
    typename blno_nnz_view_t_,
    typename bscalar_nnz_view_t_,
    typename clno_row_view_t_,
    typename clno_nnz_view_t_,
    typename cscalar_nnz_view_t_>
  void spgemm_numeric(
      KernelHandle *handle,
      typename KernelHandle::nnz_lno_t m,
      typename KernelHandle::nnz_lno_t n,
      typename KernelHandle::nnz_lno_t k,
      alno_row_view_t_ row_mapA,
      alno_nnz_view_t_ entriesA,
      ascalar_nnz_view_t_ valuesA,

      bool transposeA,
      blno_row_view_t_ row_mapB,
      blno_nnz_view_t_ entriesB,
      bscalar_nnz_view_t_ valuesB,
      bool transposeB,
      clno_row_view_t_ row_mapC,
      clno_nnz_view_t_ &entriesC,
      cscalar_nnz_view_t_ &valuesC
      ){


    typedef typename KernelHandle::SPGEMMHandleType spgemmHandleType;
    spgemmHandleType *sh = handle->get_spgemm_handle();
    if (!sh->is_symbolic_called()){
      spgemm_symbolic<KernelHandle,
                    alno_row_view_t_, alno_nnz_view_t_,
                    blno_row_view_t_, blno_nnz_view_t_,
                    clno_row_view_t_>(
          handle, m, n, k,
          row_mapA, entriesA, transposeA,
          row_mapB, entriesB, transposeB,
          row_mapC
          );
    }


    switch (sh->get_algorithm_type()){
    case SPGEMM_CUSPARSE:
      valuesC = typename cscalar_nnz_view_t_::non_const_type(Kokkos::ViewAllocateWithoutInitializing("valC"), entriesC.dimension_0());

      Impl::cuSPARSE_apply<spgemmHandleType>(
          sh,
          m,n,k,
          row_mapA, entriesA, valuesA, transposeA,
          row_mapB, entriesB, valuesB, transposeB,
          row_mapC, entriesC, valuesC);
      break;
    case SPGEMM_CUSP:
      Impl::CUSP_apply<spgemmHandleType,
        alno_row_view_t_,
        alno_nnz_view_t_,
        ascalar_nnz_view_t_,
        blno_row_view_t_,
        blno_nnz_view_t_,
        bscalar_nnz_view_t_,
        clno_row_view_t_,
        clno_nnz_view_t_,
        cscalar_nnz_view_t_ >(
          sh,
          m,n,k,
          row_mapA, entriesA, valuesA, transposeA,
          row_mapB, entriesB, valuesB, transposeB,
          row_mapC, entriesC, valuesC);
          break;
    case SPGEMM_MKL:
      Impl::mkl_apply<spgemmHandleType>(
                sh,
                m,n,k,
                row_mapA, entriesA, valuesA, transposeA,
                row_mapB, entriesB, valuesB, transposeB,
                row_mapC, entriesC, valuesC, handle->get_verbose());
      break;
    case SPGEMM_VIENNA:
      Impl::viennaCL_apply<spgemmHandleType>(
                sh,
                m,n,k,
                row_mapA, entriesA, valuesA, transposeA,
                row_mapB, entriesB, valuesB, transposeB,
                row_mapC, entriesC, valuesC, handle->get_verbose());
      break;

    case SPGEMM_DEFAULT:
    case SPGEMM_KK_MEMSPEED:
    case SPGEMM_KK_SPEED:
    case SPGEMM_KK_MEMORY:
    case SPGEMM_KK_COLOR:
    case SPGEMM_KK_MULTICOLOR:
    case SPGEMM_KK_MULTICOLOR2:
    {
      KokkosKernels::Experimental::Graph::Impl::KokkosSPGEMM
      <KernelHandle,
      alno_row_view_t_, alno_nnz_view_t_, ascalar_nnz_view_t_,
      blno_row_view_t_, blno_nnz_view_t_,  bscalar_nnz_view_t_>
      kspgemm (handle,m,n,k,row_mapA, entriesA, valuesA, transposeA, row_mapB, entriesB, valuesB, transposeB);
      kspgemm.KokkosSPGEMM_numeric(row_mapC, entriesC, valuesC);
    }
    break;

    case SPGEMM_SERIAL:
    case SPGEMM_DEBUG:
      KokkosKernels::Experimental::Graph::Impl::spgemm_debug(
          handle,
          m,
          n,
          k,
          row_mapA,
          entriesA,
          valuesA,

          transposeA,
          row_mapB,
          entriesB,
          valuesB,
          transposeB,
          row_mapC,
          entriesC,
          valuesC
          );
    default:
      break;
    }
  }
}
}
}
#endif
