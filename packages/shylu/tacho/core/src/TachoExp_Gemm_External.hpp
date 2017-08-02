#ifndef __TACHOEXP_GEMM_EXTERNAL_HPP__
#define __TACHOEXP_GEMM_EXTERNAL_HPP__


/// \file  Tacho_Gemm_External.hpp
/// \brief BLAS general matrix matrix multiplication
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Teuchos_BLAS.hpp"

namespace Tacho {

  namespace Experimental {

    template<typename ArgTransA, typename ArgTransB>
    struct Gemm<ArgTransA,ArgTransB,Algo::External> {
      template<typename PolicyType,
               typename MemberType,
               typename ScalarType,
               typename ViewTypeA,
               typename ViewTypeB,
               typename ViewTypeC>
      inline
      static int
      invoke(const PolicyType &policy,
             const MemberType &member,
             const ScalarType alpha,
             const ViewTypeA &A,
             const ViewTypeB &B,
             const ScalarType beta,
             const ViewTypeC &C) {
        
        typedef typename ViewTypeA::non_const_value_type value_type;
        typedef typename ViewTypeB::non_const_value_type value_type_b;
        typedef typename ViewTypeC::non_const_value_type value_type_c;
        
        static_assert(ViewTypeA::rank == 2,"A is not rank 2 view.");
        static_assert(ViewTypeB::rank == 2,"B is not rank 2 view.");
        static_assert(ViewTypeC::rank == 2,"C is not rank 2 view.");
        
        static_assert(std::is_same<value_type,value_type_b>::value &&
                      std::is_same<value_type_b,value_type_c>::value,
                      "A, B and C do not have the same value type.");

        const ordinal_type 
          m = C.dimension_0(),
          n = C.dimension_1(),
          k = (std::is_same<ArgTransB,Trans::NoTranspose>::value ? B.dimension_0() : B.dimension_1());
        
        if (m > 0 && n > 0 && k > 0) {
          if (get_team_rank(member) == 0) {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
#if defined( HAVE_SHYLUTACHO_MKL )
            const value_type alpha_value(alpha), beta_value(beta);
            if      (std::is_same<value_type,float>::value) 
              cblas_sgemm (CblasColMajor, ArgTransA::mkl_param, ArgTransB::mkl_param, 
                           m, n, k, 
                           (const float)alpha, 
                           (const float *)A.data(), (const MKL_INT)A.stride_1(), 
                           (const float *)B.data(), (const MKL_INT)B.stride_1(), 
                           (const float)beta, 
                           (float *)C.data(), (const MKL_INT)C.stride_1());
            else if (std::is_same<value_type,double>::value) 
              cblas_dgemm (CblasColMajor, ArgTransA::mkl_param, ArgTransB::mkl_param, 
                           m, n, k, 
                           (const double)alpha, 
                           (const double *)A.data(), (const MKL_INT)A.stride_1(), 
                           (const double *)B.data(), (const MKL_INT)B.stride_1(), 
                           (const double)beta, 
                           (double *)C.data(), (const MKL_INT)C.stride_1());
            else if (std::is_same<value_type,Kokkos::complex<float> >::value ||
                     std::is_same<value_type,   std::complex<float> >::value)
              cblas_cgemm (CblasColMajor, ArgTransA::mkl_param, ArgTransB::mkl_param, 
                           m, n, k, 
                           (const MKL_Complex8 *)&alpha_value, 
                           (const MKL_Complex8 *)A.data(), (const MKL_INT)A.stride_1(), 
                           (const MKL_Complex8 *)B.data(), (const MKL_INT)B.stride_1(), 
                           (const MKL_Complex8 *)&beta_value, 
                           (MKL_Complex8 *)C.data(), (const MKL_INT)C.stride_1());
            else if (std::is_same<value_type,Kokkos::complex<double> >::value ||
                     std::is_same<value_type,   std::complex<double> >::value)
              cblas_zgemm (CblasColMajor, ArgTransA::mkl_param, ArgTransB::mkl_param, 
                           m, n, k, 
                           (const MKL_Complex16 *)&alpha_value, 
                           (const MKL_Complex16 *)A.data(), (const MKL_INT)A.stride_1(), 
                           (const MKL_Complex16 *)B.data(), (const MKL_INT)B.stride_1(), 
                           (const MKL_Complex16 *)&beta_value, 
                           (MKL_Complex16 *)C.data(), (const MKL_INT)C.stride_1());
            else {
              TACHO_TEST_FOR_ABORT( true, ">> Datatype is not supported.");                           
            }
#else
            typedef typename TypeTraits<value_type>::std_value_type std_value_type; 
            Teuchos::BLAS<ordinal_type,std_value_type> blas;
            
            blas.GEMM(ArgTransA::teuchos_param,
                      ArgTransB::teuchos_param,
                      m, n, k,
                      std_value_type(alpha),
                      (std_value_type*)A.data(), A.stride_1(),
                      (std_value_type*)B.data(), B.stride_1(),
                      std_value_type(beta),
                      (std_value_type*)C.data(), C.stride_1());
#endif
#else
            TACHO_TEST_FOR_ABORT( true, ">> This function is only allowed in host space.");
#endif
          }
        }
        return 0;
      }
    };
  }
}
#endif
