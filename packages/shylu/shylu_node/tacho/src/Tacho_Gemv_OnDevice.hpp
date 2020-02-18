#ifndef __TACHO_GEMV_ON_DEVICE_HPP__
#define __TACHO_GEMV_ON_DEVICE_HPP__


/// \file  Tacho_Gemv_OnDevice.hpp
/// \brief BLAS general matrix matrix multiplication
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

    template<typename ArgTrans>
    struct Gemv<ArgTrans,Algo::OnDevice> {
      template<typename ScalarType,
               typename ViewTypeA,
               typename ViewTypeB,
               typename ViewTypeC>
      inline
      static int
      blas_invoke(const ScalarType alpha,
                  const ViewTypeA &A,
                  const ViewTypeB &B,
                  const ScalarType beta,
                  const ViewTypeC &C) {
        typedef typename ViewTypeA::non_const_value_type value_type;
        const ordinal_type 
          m = A.extent(0),
          n = A.extent(1),
          k = C.extent(1);

        if (m > 0 && n > 0 && k > 0) {
          for (ordinal_type p=0,offsB=0,offsC=0;p<k;++p,offsB+=B.stride_1(),offsC+=C.stride_1()) {
            Blas<value_type>::gemv(ArgTrans::param, 
                                   m, n, 
                                   value_type(alpha),
                                   A.data(), A.stride_1(),
                                   (B.data() + offsB), B.stride_0(),
                                   value_type(beta), 
                                   (C.data() + offsC), C.stride_0());
          }
        }
        return 0;
      }
#if defined(KOKKOS_ENABLE_CUDA)
      template<typename ScalarType,
               typename ViewTypeA,
               typename ViewTypeB,
               typename ViewTypeC>
      inline
      static int
      cublas_invoke(cublasHandle_t &handle,
                    const ScalarType alpha,
                    const ViewTypeA &A,
                    const ViewTypeB &B,
                    const ScalarType beta,
                    const ViewTypeC &C) {
        typedef typename ViewTypeA::non_const_value_type value_type;
        const ordinal_type 
          m = A.extent(0),
          n = A.extent(1),
          k = C.extent(1);
        int r_val(0);
        if (m > 0 && n > 0 && k > 0) {
          for (ordinal_type p=0,offsB=0,offsC=0;p<k;++p,offsB+=B.stride_1(),offsC+=C.stride_1()) {
            if (std::is_same<value_type,double>::value) 
              r_val = cublasDgemv(handle, 
                                  ArgTrans::cublas_param,
                                  m, n,
                                  &alpha,
                                  A.data(), A.stride_1(),
                                  (B.data() + offsB), B.stride_0(),
                                  &beta,
                                  (C.data() + offsC), C.stride_0());
          }
        }
        return r_val;
      }
#endif
      
      template<typename MemberType,
               typename ScalarType,
               typename ViewTypeA,
               typename ViewTypeB,
               typename ViewTypeC>
      inline
      static int
      invoke(MemberType &member,
             const ScalarType alpha,
             const ViewTypeA &A,
             const ViewTypeB &B,
             const ScalarType beta,
             const ViewTypeC &C) {
        typedef typename ViewTypeA::non_const_value_type value_type;
        typedef typename ViewTypeB::non_const_value_type value_type_b;
        typedef typename ViewTypeC::non_const_value_type value_type_c;

        typedef typename ViewTypeA::memory_space memory_space;
        typedef typename ViewTypeB::memory_space memory_space_b;
        typedef typename ViewTypeC::memory_space memory_space_c;
        
        static_assert(ViewTypeA::rank == 2,"A is not rank 2 view.");
        static_assert(ViewTypeB::rank == 2,"B is not rank 2 view.");
        static_assert(ViewTypeC::rank == 2,"C is not rank 2 view.");
        
        static_assert(std::is_same<value_type,value_type_b>::value &&
                      std::is_same<value_type_b,value_type_c>::value,
                      "A, B and C do not have the same value type.");

        static_assert(std::is_same<memory_space,memory_space_b>::value &&
                      std::is_same<memory_space_b,memory_space_c>::value,
                      "A, B and C do not have the same memory space.");
        int r_val(0);
        if (std::is_same<memory_space,Kokkos::HostSpace>::value)         
          r_val = blas_invoke(alpha, A, B, beta, C);
#if defined (KOKKOS_ENABLE_CUDA)
        if (std::is_same<memory_space,Kokkos::CudaSpace>::value ||
            std::is_same<memory_space,Kokkos::CudaUVMSpace>::value) 
          r_val = cublas_invoke(member, alpha, A, B, beta, C);
#endif
        return r_val;
      }
    };

}
#endif
