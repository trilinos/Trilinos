#ifndef __TACHOEXP_TRSV_EXTERNAL_HPP__
#define __TACHOEXP_TRSV_EXTERNAL_HPP__


/// \file  Tacho_Trsv_External.hpp
/// \brief BLAS triangular solve matrix
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Blas_External.hpp"

namespace Tacho {
  
  namespace Experimental {
    
    template<typename ArgUplo, typename ArgTransA>
    struct Trsv<ArgUplo,ArgTransA,Algo::External> {      
      template<typename PolicyType,
               typename MemberType,
               typename DiagType,
               typename ViewTypeA,
               typename ViewTypeB>
      inline
      static int
      invoke(const PolicyType &policy,
             const MemberType &member,
             const DiagType diagA,
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
            for (ordinal_type p=0,offsB=0;p<n;++p,offsB+=B.stride_1()) {  
              Blas<value_type>::trsv(ArgUplo::param, ArgTransA::param, 
                                     diagA.param, 
                                     m,
                                     A.data(), A.stride_1(), 
                                     (B.data() + offsB), B.stride_0());
            }
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
