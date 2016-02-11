#ifndef _KOKKOS_GAUSSSEIDEL_HPP
#define _KOKKOS_GAUSSSEIDEL_HPP

#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_SPGEMM_cuSPARSE_impl.hpp"
#include "KokkosKernels_SPGEMM_CUSP_impl.hpp"
#include "KokkosKernels_SPGEMM_mkl_impl.hpp"
#include "KokkosKernels_SPGEMM_impl.hpp"

namespace KokkosKernels{

namespace Experimental{

namespace Graph{

  template <typename KernelHandle,
  typename lno_row_view_t_,
  typename lno_nnz_view_t_>

  void spgemm_symbolic(
      KernelHandle *handle,
      typename KernelHandle::row_lno_t m,
      typename KernelHandle::row_lno_t n,
      typename KernelHandle::row_lno_t k,
      lno_row_view_t_ row_mapA,
      lno_nnz_view_t_ entriesA,
      bool transposeA,
      lno_row_view_t_ row_mapB,
      lno_nnz_view_t_ entriesB,
      bool transposeB,
      lno_row_view_t_ &row_mapC,
      lno_nnz_view_t_ &entriesC
      ){

    typedef typename KernelHandle::SPGEMMHandleType spgemmHandleType;
    spgemmHandleType *sh = handle->get_spgemm_handle();
    switch (sh->get_algorithm_type()){

    case SPGEMM_CUSPARSE:
      Impl::cuSPARSE_symbolic
      <spgemmHandleType,
      lno_row_view_t_,
      lno_nnz_view_t_>(sh, m,n,k,
          row_mapA, entriesA, transposeA,
          row_mapB, entriesB, transposeB,
          row_mapC, entriesC);
      break;

    case SPGEMM_CUSP:
      break;

    case SPGEMM_KK1:
    {
      KokkosKernels::Experimental::Graph::Impl::KokkosSPGEMM
      <KernelHandle,
      lno_row_view_t_, lno_nnz_view_t_, typename KernelHandle::in_scalar_nnz_view_t,
      lno_row_view_t_, lno_nnz_view_t_, typename KernelHandle::in_scalar_nnz_view_t>
      kspgemm (handle,m,n,k,row_mapA, entriesA, transposeA, row_mapB, entriesB, transposeB);
      kspgemm.KokkosSPGEMM_symbolic(row_mapC, entriesC);
    }
      break;

    case SPGEMM_DEFAULT:
    case SPGEMM_SERIAL:
    case SPGEMM_MKL:
    default:
      break;
    }
    sh->set_call_symbolic();

  }


  template <typename KernelHandle,
    typename lno_row_view_t_,
    typename lno_nnz_view_t_,
    typename scalar_nnz_view_t_>
  void spgemm_numeric(

      KernelHandle *handle,
      typename KernelHandle::row_lno_t m,
      typename KernelHandle::row_lno_t n,
      typename KernelHandle::row_lno_t k,
      lno_row_view_t_ row_mapA,
      lno_nnz_view_t_ entriesA,
      scalar_nnz_view_t_ valuesA,

      bool transposeA,
      lno_row_view_t_ row_mapB,
      lno_nnz_view_t_ entriesB,
      scalar_nnz_view_t_ valuesB,
      bool transposeB,
      lno_row_view_t_ &row_mapC,
      lno_nnz_view_t_ &entriesC,
      scalar_nnz_view_t_ &valuesC
      ){


    typedef typename KernelHandle::SPGEMMHandleType spgemmHandleType;
    spgemmHandleType *sh = handle->get_spgemm_handle();
    if (!sh->is_symbolic_called()){
      spgemm_symbolic<KernelHandle,lno_row_view_t_, lno_nnz_view_t_>(
          handle, m, n, k,
          row_mapA, entriesA, transposeA,
          row_mapB, entriesB, transposeB,
          row_mapC, entriesC
          );
    }

    valuesC = typename scalar_nnz_view_t_::non_const_type("valC", entriesC.dimension_0());
    switch (sh->get_algorithm_type()){
    case SPGEMM_CUSPARSE:
        //no op.
      break;
    case SPGEMM_CUSP:
        //no op.
          break;

    case SPGEMM_KK1:
    case SPGEMM_DEFAULT:
    case SPGEMM_SERIAL:
    case SPGEMM_MKL:
    default:
      break;
    }

    sh->set_call_numeric();
  }

  template <typename KernelHandle,
    typename lno_row_view_t_,
    typename lno_nnz_view_t_,
    typename scalar_nnz_view_t_>
  void spgemm_apply(
      KernelHandle *handle,
      typename KernelHandle::row_lno_t m,
      typename KernelHandle::row_lno_t n,
      typename KernelHandle::row_lno_t k,
      lno_row_view_t_ row_mapA,
      lno_nnz_view_t_ entriesA,
      scalar_nnz_view_t_ valuesA,

      bool transposeA,
      lno_row_view_t_ row_mapB,
      lno_nnz_view_t_ entriesB,
      scalar_nnz_view_t_ valuesB,
      bool transposeB,
      lno_row_view_t_ &row_mapC,
      lno_nnz_view_t_ &entriesC,
      scalar_nnz_view_t_ &valuesC
      ){


    typedef typename KernelHandle::SPGEMMHandleType spgemmHandleType;
    spgemmHandleType *sh = handle->get_spgemm_handle();
    if (!sh->is_numeric_called()){
      spgemm_numeric<KernelHandle,
        lno_row_view_t_,
        lno_nnz_view_t_,
        scalar_nnz_view_t_>(
          handle, m, n, k,
          row_mapA, entriesA, valuesA, transposeA,
          row_mapB, entriesB, valuesB, transposeB,
          row_mapC, entriesC, valuesC
          );
    }

    //valuesC = typename KernelHandle::value_array_type("valC", entriesC.dimension_0());
    switch (sh->get_algorithm_type()){
    case SPGEMM_CUSPARSE:
      std::cout << "SPGEMM_CUSPARSE" << std::endl;
      Impl::cuSPARSE_apply<spgemmHandleType>(
          sh,
          m,n,k,
          row_mapA, entriesA, valuesA, transposeA,
          row_mapB, entriesB, valuesB, transposeB,
          row_mapC, entriesC, valuesC);
      break;
    case SPGEMM_CUSP:
      std::cout << "SPGEMM_CUSP" << std::endl;
      Impl::CUSP_apply<spgemmHandleType,
        lno_row_view_t_,
        lno_nnz_view_t_,
        scalar_nnz_view_t_>(
          sh,
          m,n,k,
          row_mapA, entriesA, valuesA, transposeA,
          row_mapB, entriesB, valuesB, transposeB,
          row_mapC, entriesC, valuesC);
          break;
    case SPGEMM_MKL:
      std::cout << "MKL" << std::endl;
      Impl::mkl_apply<spgemmHandleType>(
                sh,
                m,n,k,
                row_mapA, entriesA, valuesA, transposeA,
                row_mapB, entriesB, valuesB, transposeB,
                row_mapC, entriesC, valuesC);
      break;

    case SPGEMM_KK1:
    case SPGEMM_DEFAULT:
    case SPGEMM_SERIAL:
    default:
      break;
    }
  }
}
}
}
#endif
