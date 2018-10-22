#ifndef __TACHO_GEMV_INTERNAL_HPP__
#define __TACHO_GEMV_INTERNAL_HPP__


/// \file  Tacho_Gemv_Internal.hpp
/// \brief BLAS general matrix matrix multiplication
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Blas_Team.hpp"

namespace Tacho {

    template<typename ArgTrans>
    struct Gemv<ArgTrans,Algo::Internal> {
      template<typename PolicyType,
               typename MemberType,
               typename ScalarType,
               typename ViewTypeA,
               typename ViewTypeB,
               typename ViewTypeC>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(PolicyType &policy,
             MemberType &member,
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

        const ordinal_type m = A.extent(0);
        const ordinal_type n = A.extent(1);
        const ordinal_type k = C.extent(1);

        if (m > 0 && n > 0 && k > 0) 
          for (ordinal_type p=0,offsB=0,offsC=0;p<k;++p,offsB+=B.stride_1(),offsC+=C.stride_1()) {
            BlasTeam<value_type>::gemv(member,
                                       ArgTrans::param, 
                                       m, n, 
                                       value_type(alpha),
                                       A.data(), A.stride_1(),
                                       (B.data() + offsB), B.stride_0(),
                                       value_type(beta), 
                                       (C.data() + offsC), C.stride_0());
          }
        return 0;
      }
    };

}
#endif
