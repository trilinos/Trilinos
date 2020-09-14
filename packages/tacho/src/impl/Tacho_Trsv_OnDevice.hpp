#ifndef __TACHO_TRSV_ON_DEVICE_HPP__
#define __TACHO_TRSV_ON_DEVICE_HPP__


/// \file  Tacho_Trsv_OnDevice.hpp
/// \brief BLAS triangular solve matrix
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Blas_External.hpp"

namespace Tacho {
  
    template<typename ArgUplo, typename ArgTransA>
    struct Trsv<ArgUplo,ArgTransA,Algo::OnDevice> {      
      template<typename DiagType,
               typename ViewTypeA,
               typename ViewTypeB>
      inline
      static int
      blas_invoke(const DiagType diagA,
                  const ViewTypeA &A,
                  const ViewTypeB &B) {
        typedef typename ViewTypeA::non_const_value_type value_type;        
        const ordinal_type m = B.extent(0), n = B.extent(1);
        if (m > 0 && n > 0) {
          if (n == 1) {
            Blas<value_type>::trsv(ArgUplo::param, ArgTransA::param, 
                                   diagA.param, 
                                   m,
                                   A.data(), A.stride_1(), 
                                   B.data(), B.stride_0());
          } else {
            Blas<value_type>::trsm(Side::Left::param, 
                                   ArgUplo::param, 
                                   ArgTransA::param, 
                                   diagA.param,
                                   m, n,
                                   value_type(1),
                                   A.data(), A.stride_1(),
                                   B.data(), B.stride_1());
          }
        }
        return 0;
      }

#if defined (KOKKOS_ENABLE_CUDA)
      template<typename DiagType,
               typename ViewTypeA,
               typename ViewTypeB>
      inline
      static int
      cublas_invoke(cublasHandle_t &handle,
                    const DiagType diagA,
                    const ViewTypeA &A,
                    const ViewTypeB &B) {
        typedef typename ViewTypeA::non_const_value_type value_type;        
        const ordinal_type m = B.extent(0), n = B.extent(1);
        int r_val(0); 
        if (m > 0 && n > 0) {
          if (n == 1) {
            r_val = Blas<value_type>::trsv(handle,
                                           ArgUplo::cublas_param,
                                           ArgTransA::cublas_param,
                                           diagA.cublas_param,
                                           m, 
                                           A.data(), A.stride_1(),
                                           B.data(), B.stride_0()); 
          } else {
            r_val = Blas<value_type>::trsm(handle,
                                           Side::Left::cublas_param,
                                           ArgUplo::cublas_param,
                                           ArgTransA::cublas_param,
                                           diagA.cublas_param,
                                           m, n, 
                                           value_type(1),
                                           A.data(), A.stride_1(),
                                           B.data(), B.stride_1()); 
          }
        }
        return r_val;
      }
#endif

      template<typename MemberType,
               typename DiagType,
               typename ViewTypeA,
               typename ViewTypeB>
      inline
      static int
      invoke(MemberType &member,
             const DiagType diagA,
             const ViewTypeA &A,
             const ViewTypeB &B) {
        typedef typename ViewTypeA::non_const_value_type value_type;
        typedef typename ViewTypeB::non_const_value_type value_type_b;

        typedef typename ViewTypeA::memory_space memory_space;
        typedef typename ViewTypeB::memory_space memory_space_b;
        
        static_assert(ViewTypeA::rank == 2,"A is not rank 2 view.");
        static_assert(ViewTypeB::rank == 2,"B is not rank 2 view.");
        
        static_assert(std::is_same<value_type,value_type_b>::value,
                      "A and B do not have the same value type.");

        static_assert(std::is_same<memory_space,memory_space_b>::value,
                      "A and B do not have the same memory space.");
        int r_val(0);
        if (std::is_same<memory_space,Kokkos::HostSpace>::value) 
          r_val = blas_invoke(diagA, A, B);
#if defined (KOKKOS_ENABLE_CUDA)
        if (std::is_same<memory_space,Kokkos::CudaSpace>::value ||
            std::is_same<memory_space,Kokkos::CudaUVMSpace>::value) 
          r_val = cublas_invoke(member, diagA, A, B);
#endif
        return r_val;
      }
    };

}
#endif
