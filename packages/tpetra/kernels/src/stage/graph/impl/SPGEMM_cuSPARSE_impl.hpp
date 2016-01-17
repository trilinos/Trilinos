
#ifndef _KOKKOSSPGEMMCUSPARSE_HPP
#define _KOKKOSSPGEMMCUSPARSE_HPP

//#define KERNELS_HAVE_CUSPARSE

#ifdef KERNELS_HAVE_CUSPARSE
#include "cusparse.h"
#endif
namespace KokkosKernels{

namespace Experimental{

namespace Graph{
namespace Impl{

  template <typename KernelHandle>
  void cuSPARSE_symbolic(
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
      typename KernelHandle::idx_array_type &entriesC){

#ifdef KERNELS_HAVE_CUSPARSE

    typedef typename KernelHandle::idx_device_type device1;
    typedef typename KernelHandle::idx_edge_device_type device2;
    typedef typename KernelHandle::idx idx;
    typedef typename KernelHandle::idx_array_type idx_array_type;


    if (Kokkos::Impl::is_same<Kokkos::Cuda, device1 >::value){
      std::cerr << "MEMORY IS NOT ALLOCATED IN GPU DEVICE for CUSPARSE" << std::endl;
      return;
    }
    if (Kokkos::Impl::is_same<Kokkos::Cuda, device2 >::value){
      std::cerr << "MEMORY IS NOT ALLOCATED IN GPU DEVICE for CUSPARSE" << std::endl;
      return;
    }

    if (Kokkos::Impl::is_same<idx, int>::value){
      row_mapC = idx_array_type("rowMapC", m + 1);
      idx *a_xadj = row_mapA.ptr_on_device();
      idx *b_xadj = row_mapB.ptr_on_device();
      idx *c_xadj = row_mapC.ptr_on_device();

      idx *a_adj = entriesA.ptr_on_device();
      idx *b_adj = entriesB.ptr_on_device();
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
      entriesC = idx_array_type("entriesC", nnzC);
    }
    else {
      std::cerr << "CUSPARSE requires integer values" << std::endl;
      return;
    }
#else
    std::cerr << "CUSPARSE IS NOT DEFINED" << std::endl;
    return;
#endif

  }


  template <typename KernelHandle>
  void cuSPARSE_apply(
      KernelHandle *handle,
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
      typename KernelHandle::value_array_type &valuesC){

#ifdef KERNELS_HAVE_CUSPARSE
    typedef typename KernelHandle::idx idx;
    typedef typename KernelHandle::idx_array_type idx_array_type;

    typedef typename KernelHandle::value_type value_type;


    typedef typename KernelHandle::idx_device_type device1;
    typedef typename KernelHandle::idx_edge_device_type device2;
    typedef typename KernelHandle::value_type_device_type device3;

    std::cout << "RUNNING CUSParse" << std::endl;

    if (Kokkos::Impl::is_same<Kokkos::Cuda, device1 >::value){
      std::cerr << "MEMORY IS NOT ALLOCATED IN GPU DEVICE for CUSPARSE" << std::endl;
      return;
    }
    if (Kokkos::Impl::is_same<Kokkos::Cuda, device2 >::value){
      std::cerr << "MEMORY IS NOT ALLOCATED IN GPU DEVICE for CUSPARSE" << std::endl;
      return;
    }
    if (Kokkos::Impl::is_same<Kokkos::Cuda, device3 >::value){
      std::cerr << "MEMORY IS NOT ALLOCATED IN GPU DEVICE for CUSPARSE" << std::endl;
      return;
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

      value_type *a_ew = valuesA.ptr_on_device();
      value_type *b_ew = valuesB.ptr_on_device();
      value_type *c_ew = valuesC.ptr_on_device();

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
        std::cerr << "CUSPARSE requires float or double values. cuComplex and cuDoubleComplex are not implemented yet." << std::endl;
        return;
      }




    }
    else {
      std::cerr << "CUSPARSE requires integer values" << std::endl;
      return;
    }
#else
    std::cerr << "CUSPARSE IS NOT DEFINED" << std::endl;
    return;
#endif
  }
}
}
}
}

#endif
