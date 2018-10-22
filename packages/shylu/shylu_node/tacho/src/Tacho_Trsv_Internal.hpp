#ifndef __TACHO_TRSV_INTERNAL_HPP__
#define __TACHO_TRSV_INTERNAL_HPP__


/// \file  Tacho_Trsv_Internal.hpp
/// \brief BLAS triangular solve matrix
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Blas_Team.hpp"

namespace Tacho {
  
    template<typename ArgUplo, typename ArgTransA>
    struct Trsv<ArgUplo,ArgTransA,Algo::Internal> {      
      template<typename PolicyType,
               typename MemberType,
               typename DiagType,
               typename ViewTypeA,
               typename ViewTypeB>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(PolicyType &policy,
             MemberType &member,
             const DiagType diagA,
             const ViewTypeA &A,
             const ViewTypeB &B) {
        typedef typename ViewTypeA::non_const_value_type value_type;
        typedef typename ViewTypeB::non_const_value_type value_type_b;
        
        static_assert(ViewTypeA::rank == 2,"A is not rank 2 view.");
        static_assert(ViewTypeB::rank == 2,"B is not rank 2 view.");
        
        static_assert(std::is_same<value_type,value_type_b>::value,
                      "A and B do not have the same value type.");
        
        const ordinal_type m = B.extent(0);
        const ordinal_type n = B.extent(1);
        
        if (m > 0 && n > 0) 
          for (ordinal_type p=0,offsB=0;p<n;++p,offsB+=B.stride_1()) {  
            BlasTeam<value_type>::trsv(member,
                                       ArgUplo::param, ArgTransA::param, 
                                       diagA.param, 
                                       m,
                                       A.data(), A.stride_1(), 
                                       (B.data() + offsB), B.stride_0());
          }
        return 0;
      }
    };

}
#endif
