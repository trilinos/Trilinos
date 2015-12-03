#ifndef _KOKKOS_GAUSSSEIDEL_HPP
#define _KOKKOS_GAUSSSEIDEL_HPP

#include "KokkosKernelsHandle.hpp"
#include "SPGEMM_cuSPARSE_impl.hpp"
#include "SPGEMM_CUSP_impl.hpp"
namespace KokkosKernels{

namespace Experimental{

namespace Graph{

  template <class KernelHandle>
  void spgemm_symbolic(
      KernelHandle *handle,
      typename KernelHandle::idx m,
      typename KernelHandle::idx n,
      typename KernelHandle::idx k,
      typename KernelHandle::idx_array_type row_mapA,
      typename KernelHandle::idx_edge_array_type entriesA,
      bool transposeA,
      typename KernelHandle::idx_array_type row_mapB,
      typename KernelHandle::idx_edge_array_type entriesB,
      bool transposeB,
      typename KernelHandle::idx_array_type &row_mapC,
      typename KernelHandle::idx_array_type &entriesC
      ){

    typedef typename KernelHandle::SPGEMMHandleType spgemmHandleType;
    spgemmHandleType *sh = handle->get_spgemm_handle();
    switch (sh->get_algorithm_type()){
    case SPGEMM_CUSPARSE:
      Impl::cuSPARSE_symbolic<spgemmHandleType>(sh, m,n,k,
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

  template <class KernelHandle>
  void spgemm_numeric(KernelHandle *handle,
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
      typename KernelHandle::value_array_type &valuesC
      ){


    typedef typename KernelHandle::SPGEMMHandleType spgemmHandleType;
    spgemmHandleType *sh = handle->get_spgemm_handle();
    if (!sh->is_symbolic_called()){
      spgemm_symbolic<KernelHandle>(
          handle, m, n, k,
          row_mapA, entriesA, transposeA,
          row_mapB, entriesB, transposeB,
          row_mapC, entriesC
          );
    }

    valuesC = typename KernelHandle::value_array_type("valC", entriesC.dimension_0());
    switch (sh->get_algorithm_type()){
    case SPGEMM_CUSPARSE:
      Impl::cuSPARSE_numeric<spgemmHandleType>(
          sh,
          m,n,k,
          row_mapA, entriesA, valuesA, transposeA,
          row_mapB, entriesB, valuesB, transposeB,
          row_mapC, entriesC, valuesC);
      break;
    case SPGEMM_CUSP:
      Impl::CUSP_numeric<spgemmHandleType>(
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
