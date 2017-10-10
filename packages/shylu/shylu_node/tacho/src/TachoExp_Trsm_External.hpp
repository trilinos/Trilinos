#ifndef __TACHOEXP_TRSM_EXTERNAL_HPP__
#define __TACHOEXP_TRSM_EXTERNAL_HPP__


/// \file  Tacho_Trsm_External.hpp
/// \brief BLAS triangular solve matrix
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Blas_External.hpp"

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
#if defined( TACHO_PROFILE_TIME_PER_THREAD )
            Kokkos::Impl::Timer timer;
            timer.reset();
#endif
            Blas<value_type>::trsm(ArgSide::param, 
                                   ArgUplo::param, 
                                   ArgTransA::param, 
                                   diagA.param,
                                   m, n,
                                   value_type(alpha),
                                   A.data(), A.stride_1(),
                                   B.data(), B.stride_1());
#if defined( TACHO_PROFILE_TIME_PER_THREAD )
            g_time_per_thread[omp_get_thread_num()] += timer.seconds();
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
