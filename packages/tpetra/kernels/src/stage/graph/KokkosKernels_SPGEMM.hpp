#ifndef _KOKKOS_GAUSSSEIDEL_HPP
#define _KOKKOS_GAUSSSEIDEL_HPP

#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_SPGEMM_cuSPARSE_impl.hpp"
#include "KokkosKernels_SPGEMM_CUSP_impl.hpp"
#include "KokkosKernels_SPGEMM_mkl_impl.hpp"

namespace KokkosKernels{

namespace Experimental{

namespace Graph{

  template <typename KernelHandle,
  typename in_row_index_view_type,
  typename in_nonzero_index_view_type>

  void spgemm_symbolic(
      KernelHandle *handle,
      typename KernelHandle::row_index_type m,
      typename KernelHandle::row_index_type n,
      typename KernelHandle::row_index_type k,
      in_row_index_view_type row_mapA,
      in_nonzero_index_view_type entriesA,
      bool transposeA,
      in_row_index_view_type row_mapB,
      in_nonzero_index_view_type entriesB,
      bool transposeB,
      in_row_index_view_type &row_mapC,
      in_nonzero_index_view_type &entriesC
      ){

    typedef typename KernelHandle::SPGEMMHandleType spgemmHandleType;
    spgemmHandleType *sh = handle->get_spgemm_handle();
    switch (sh->get_algorithm_type()){
    case SPGEMM_CUSPARSE:
      Impl::cuSPARSE_symbolic
        <spgemmHandleType,
          in_row_index_view_type,
          in_nonzero_index_view_type>(sh, m,n,k,
          row_mapA, entriesA, transposeA,
          row_mapB, entriesB, transposeB,
          row_mapC, entriesC);
      break;
    case SPGEMM_CUSP:
      break;
    default:
      break;
    }
    sh->set_call_symbolic();

  }


  template <typename KernelHandle,
    typename in_row_index_view_type,
    typename in_nonzero_index_view_type,
    typename in_nonzero_value_view_type>
  void spgemm_numeric(

      KernelHandle *handle,
      typename KernelHandle::row_index_type m,
      typename KernelHandle::row_index_type n,
      typename KernelHandle::row_index_type k,
      in_row_index_view_type row_mapA,
      in_nonzero_index_view_type entriesA,
      in_nonzero_value_view_type valuesA,

      bool transposeA,
      in_row_index_view_type row_mapB,
      in_nonzero_index_view_type entriesB,
      in_nonzero_value_view_type valuesB,
      bool transposeB,
      in_row_index_view_type &row_mapC,
      in_nonzero_index_view_type &entriesC,
      in_nonzero_value_view_type &valuesC
      ){


    typedef typename KernelHandle::SPGEMMHandleType spgemmHandleType;
    spgemmHandleType *sh = handle->get_spgemm_handle();
    if (!sh->is_symbolic_called()){
      spgemm_symbolic<KernelHandle,in_row_index_view_type, in_nonzero_index_view_type>(
          handle, m, n, k,
          row_mapA, entriesA, transposeA,
          row_mapB, entriesB, transposeB,
          row_mapC, entriesC
          );
    }

    valuesC = typename in_nonzero_value_view_type::non_const_type("valC", entriesC.dimension_0());
    switch (sh->get_algorithm_type()){
    case SPGEMM_CUSPARSE:
        //no op.
      break;
    case SPGEMM_CUSP:
        //no op.
          break;
    default:
      break;
    }

    sh->set_call_numeric();
  }

  template <typename KernelHandle,
    typename in_row_index_view_type,
    typename in_nonzero_index_view_type,
    typename in_nonzero_value_view_type>
  void spgemm_apply(
      KernelHandle *handle,
      typename KernelHandle::row_index_type m,
      typename KernelHandle::row_index_type n,
      typename KernelHandle::row_index_type k,
      in_row_index_view_type row_mapA,
      in_nonzero_index_view_type entriesA,
      in_nonzero_value_view_type valuesA,

      bool transposeA,
      in_row_index_view_type row_mapB,
      in_nonzero_index_view_type entriesB,
      in_nonzero_value_view_type valuesB,
      bool transposeB,
      in_row_index_view_type &row_mapC,
      in_nonzero_index_view_type &entriesC,
      in_nonzero_value_view_type &valuesC
      ){


    typedef typename KernelHandle::SPGEMMHandleType spgemmHandleType;
    spgemmHandleType *sh = handle->get_spgemm_handle();
    if (!sh->is_numeric_called()){
      spgemm_numeric<KernelHandle,
        in_row_index_view_type,
        in_nonzero_index_view_type,
        in_nonzero_value_view_type>(
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
        in_row_index_view_type,
        in_nonzero_index_view_type,
        in_nonzero_value_view_type>(
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
    default:
      break;
    }
  }
}
}
}
#endif
