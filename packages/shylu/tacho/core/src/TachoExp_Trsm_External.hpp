#ifndef __TACHOEXP_TRSM_EXTERNAL_HPP__
#define __TACHOEXP_TRSM_EXTERNAL_HPP__


/// \file  Tacho_Trsm_External.hpp
/// \brief BLAS triangular solve matrix
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Teuchos_BLAS.hpp"

namespace Tacho {
  
  namespace Experimental {
    
    template<typename ArgSide, typename ArgUplo, typename ArgTransA>
    struct Trsm<ArgSide,ArgUplo,ArgTransA,Algo::External> {      
      template<typename PolicyType,
               typename MemberType,
               typename DiagType,
               typename ScalarType,
               typename ViewTypeA,
               typename ViewTypeB>
      inline
      static int
      invoke(const PolicyType &policy,
             const MemberType &member,
             const DiagType diagA,
             const ScalarType alpha,
             const ViewTypeA &A,
             const ViewTypeB &B) {
        typedef typename ViewTypeA::non_const_value_type value_type;
        typedef typename ViewTypeB::non_const_value_type value_type_b;
        
        static_assert(ViewTypeA::rank == 2,"A is not rank 2 view.");
        static_assert(ViewTypeB::rank == 2,"B is not rank 2 view.");
        
        static_assert(std::is_same<value_type,value_type_b>::value,
                      "A and B do not have the same value type.");
        
        const ordinal_type m = B.dimension_0(), n = B.dimension_1();
        
        if (m > 0 && n > 0) {
          if (get_team_rank(member) == 0) {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
#if defined( HAVE_SHYLUTACHO_MKL )
            const value_type alpha_value(alpha);
            if      (std::is_same<value_type,float>::value) 
              cblas_strsm (CblasColMajor, ArgSide::mkl_param, ArgUplo::mkl_param, ArgTransA::mkl_param, 
                           diagA.mkl_param, 
                           m, n, 
                           (const float)alpha, 
                           (const float *)A.data(), (const MKL_INT)A.stride_1(), 
                           (float *)B.data(), (const MKL_INT)B.stride_1());
            else if (std::is_same<value_type,double>::value) 
              cblas_dtrsm (CblasColMajor, ArgSide::mkl_param, ArgUplo::mkl_param, ArgTransA::mkl_param, 
                           diagA.mkl_param, 
                           m, n, 
                           (const double)alpha, 
                           (const double *)A.data(), (const MKL_INT)A.stride_1(), 
                           (double *)B.data(), (const MKL_INT)B.stride_1());
            else if (std::is_same<value_type,Kokkos::complex<float> >::value ||
                     std::is_same<value_type,   std::complex<float> >::value)
              cblas_ctrsm (CblasColMajor, ArgSide::mkl_param, ArgUplo::mkl_param, ArgTransA::mkl_param, 
                           diagA.mkl_param, 
                           m, n, 
                           (const MKL_Complex8 *)&alpha_value, 
                           (const MKL_Complex8 *)A.data(), (const MKL_INT)A.stride_1(), 
                           (MKL_Complex8 *)B.data(), (const MKL_INT)B.stride_1());
            else if (std::is_same<value_type,Kokkos::complex<double> >::value ||
                     std::is_same<value_type,   std::complex<double> >::value)
              cblas_ztrsm (CblasColMajor, ArgSide::mkl_param, ArgUplo::mkl_param, ArgTransA::mkl_param, 
                           diagA.mkl_param, 
                           m, n, 
                           (const MKL_Complex16 *)&alpha_value, 
                           (const MKL_Complex16 *)A.data(), (const MKL_INT)A.stride_1(), 
                           (MKL_Complex16 *)B.data(), (const MKL_INT)B.stride_1());
            else {
              TACHO_TEST_FOR_ABORT( true, ">> Datatype is not supported.");                           
            }
#else
            typedef typename TypeTraits<value_type>::std_value_type std_value_type;
            Teuchos::BLAS<ordinal_type,std_value_type> blas;
            blas.TRSM(ArgSide::teuchos_param, 
                      ArgUplo::teuchos_param, 
                      ArgTransA::teuchos_param, 
                      diagA.teuchos_param,
                      m, n,
                      std_value_type(alpha),
                      (std_value_type*)A.data(), A.stride_1(),
                      (std_value_type*)B.data(), B.stride_1());
#endif
#else
            TACHO_TEST_FOR_ABORT( true, "This function is only allowed in host space.");
#endif
          }
        }
        return 0;
      }
    };
  }
}
#endif
