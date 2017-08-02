#ifndef __TACHOEXP_HERK_EXTERNAL_HPP__
#define __TACHOEXP_HERK_EXTERNAL_HPP__


/// \file  Tacho_Herk_External.hpp
/// \brief BLAS hermitian rank-k update
/// \author Kyungjoo Kim (kyukim@sandia.gov)
#include "Teuchos_BLAS.hpp"

namespace Tacho {

  namespace Experimental {

    template<typename ArgUplo, typename ArgTrans>
    struct Herk<ArgUplo,ArgTrans,Algo::External> {
      template<typename SchedType,
               typename MemberType,
               typename ScalarType,
               typename ViewTypeA,
               typename ViewTypeC>
      inline
      static int
      invoke(const SchedType &sched,
             const MemberType &member,
             const ScalarType alpha,
             const ViewTypeA &A,
             const ScalarType beta,
             const ViewTypeC &C) {
        
        typedef typename ViewTypeA::non_const_value_type value_type;
        typedef typename ViewTypeC::non_const_value_type value_type_c;
        
        static_assert(ViewTypeA::rank == 2,"A is not rank 2 view.");
        static_assert(ViewTypeC::rank == 2,"B is not rank 2 view.");
        
        static_assert(std::is_same<value_type,value_type_c>::value,
                      "A and C do not have the same value type.");
        
        const ordinal_type 
          n = C.dimension_0(), 
          k = (std::is_same<ArgTrans,Trans::NoTranspose>::value ? A.dimension_1() : A.dimension_0());
        if (n > 0 && k > 0) {
          if (get_team_rank(member) == 0) {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
#if defined( HAVE_SHYLUTACHO_MKL )
            if      (std::is_same<value_type,float>::value) 
              cblas_ssyrk (CblasColMajor, ArgUplo::mkl_param, ArgTrans::mkl_param, 
                           n, k, 
                           (const float)alpha, 
                           (const float *)A.data(), (const MKL_INT)A.stride_1(), 
                           (const float)beta, 
                           (float *)C.data(), (const MKL_INT)C.stride_1());
            else if (std::is_same<value_type,double>::value) 
              cblas_dsyrk (CblasColMajor, ArgUplo::mkl_param, ArgTrans::mkl_param, 
                           n, k, 
                           (const double)alpha, 
                           (const double *)A.data(), (const MKL_INT)A.stride_1(), 
                           (const double)beta, 
                           (double *)C.data(), (const MKL_INT)C.stride_1());
            else if (std::is_same<value_type,Kokkos::complex<float> >::value ||
                     std::is_same<value_type,   std::complex<float> >::value)
              cblas_cherk (CblasColMajor, ArgUplo::mkl_param, ArgTrans::mkl_param, 
                           n, k, 
                           (const float)real(alpha), 
                           (const MKL_Complex8 *)A.data(), (const MKL_INT)A.stride_1(), 
                           (const float)real(beta), 
                           (MKL_Complex8 *)C.data(), (const MKL_INT)C.stride_1());              
            else if (std::is_same<value_type,Kokkos::complex<double> >::value ||
                     std::is_same<value_type,   std::complex<double> >::value)
              cblas_zherk (CblasColMajor, ArgUplo::mkl_param, ArgTrans::mkl_param, 
                           n, k, 
                           (const double)real(alpha), 
                           (const MKL_Complex16 *)A.data(), (const MKL_INT)A.stride_1(), 
                           (const double)real(beta), 
                           (MKL_Complex16 *)C.data(), (const MKL_INT)C.stride_1());              
            else {
              TACHO_TEST_FOR_ABORT( true, ">> Datatype is not supported.");                           
            }
#else
            typedef typename TypeTraits<value_type>::std_value_type std_value_type; 
            Teuchos::BLAS<ordinal_type,std_value_type> blas;
            
            blas.HERK(ArgUplo::teuchos_param,
                      ArgTrans::teuchos_param,
                      n, k,
                      std_value_type(alpha),
                      (std_value_type*)A.data(), A.stride_1(),
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
