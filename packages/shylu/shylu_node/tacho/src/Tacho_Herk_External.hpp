#ifndef __TACHO_HERK_EXTERNAL_HPP__
#define __TACHO_HERK_EXTERNAL_HPP__


/// \file  Tacho_Herk_External.hpp
/// \brief BLAS hermitian rank-k update
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Blas_External.hpp"

namespace Tacho {

    template<typename ArgUplo, typename ArgTrans>
    struct Herk<ArgUplo,ArgTrans,Algo::External> {
      template<typename ScalarType,
               typename ViewTypeA,
               typename ViewTypeC>
      inline
      static int
      invoke(const ScalarType alpha,
             const ViewTypeA &A,
             const ScalarType beta,
             const ViewTypeC &C) {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )        
        typedef typename ViewTypeA::non_const_value_type value_type;
        typedef typename ViewTypeC::non_const_value_type value_type_c;
        
        static_assert(ViewTypeA::rank == 2,"A is not rank 2 view.");
        static_assert(ViewTypeC::rank == 2,"B is not rank 2 view.");
        
        static_assert(std::is_same<value_type,value_type_c>::value,
                      "A and C do not have the same value type.");
        
        const ordinal_type 
          n = C.extent(0), 
          k = (std::is_same<ArgTrans,Trans::NoTranspose>::value ? A.extent(1) : A.extent(0));
        if (n > 0 && k > 0) {
          Blas<value_type>::herk(ArgUplo::param,
                                 ArgTrans::param,
                                 n, k,
                                 value_type(alpha),
                                 A.data(), A.stride_1(),
                                 value_type(beta),
                                 C.data(), C.stride_1());
        }
#else
        TACHO_TEST_FOR_ABORT( true, ">> This function is only allowed in host space.");
#endif
        return 0;
      }

      template<typename SchedulerType,
               typename MemberType,
               typename ScalarType,
               typename ViewTypeA,
               typename ViewTypeC>
      inline
      static int
      invoke(SchedulerType &sched,
             MemberType &member,
             const ScalarType alpha,
             const ViewTypeA &A,
             const ScalarType beta,
             const ViewTypeC &C) {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )        
        //Kokkos::single(Kokkos::PerTeam(member), [&]() {
        invoke(alpha, A, beta, C);
        //});
#else
        TACHO_TEST_FOR_ABORT( true, ">> This function is only allowed in host space.");
#endif
        return 0;
      }
    };

}
#endif
