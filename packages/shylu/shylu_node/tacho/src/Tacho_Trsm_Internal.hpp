#ifndef __TACHO_TRSM_INTERNAL_HPP__
#define __TACHO_TRSM_INTERNAL_HPP__


/// \file  Tacho_Trsm_Internal.hpp
/// \brief BLAS triangular solve matrix
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Blas_Team.hpp"

namespace Tacho {
  
    template<typename ArgSide, typename ArgUplo, typename ArgTransA>
    struct Trsm<ArgSide,ArgUplo,ArgTransA,Algo::Internal> {      
      template<typename PolicyType,
               typename MemberType,
               typename DiagType,
               typename ScalarType,
               typename ViewTypeA,
               typename ViewTypeB>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(PolicyType &policy,
             MemberType &member,
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
        
        const ordinal_type m = B.dimension_0();
        const ordinal_type n = B.dimension_1();

        if (m > 0 && n > 0) 
          BlasTeam<value_type>::trsm(member,
                                     ArgSide::param, 
                                     ArgUplo::param, 
                                     ArgTransA::param, 
                                     diagA.param,
                                     m, n,
                                     value_type(alpha),
                                     A.data(), A.stride_1(),
                                     B.data(), B.stride_1());
        return 0;
      }
    };

}
#endif
