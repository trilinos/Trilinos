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
      template<typename PolicyType,
               typename MemberType,
               typename ScalarType,
               typename ViewTypeA,
               typename ViewTypeC>
      inline
      static int
      invoke(const PolicyType &policy,
             const MemberType &member,
             const ScalarType alpha,
             const ViewTypeA &A,
             const ScalarType beta,
             const ViewTypeC &C) {
        
        typedef typename ViewTypeA::non_const_value_type value_type;
        typedef typename ViewTypeC::non_const_value_type value_type_c;
        typedef typename ViewTypeA::array_layout array_layout_a;
        typedef typename ViewTypeC::array_layout array_layout_c;
        
        static_assert(ViewTypeA::rank == 2,"A is not rank 2 view.");
        static_assert(ViewTypeC::rank == 2,"B is not rank 2 view.");
        
        static_assert(std::is_same<array_layout_a,Kokkos::LayoutLeft>::value,
                      "A does not have Kokkos::LayoutLeft.");
        static_assert(std::is_same<array_layout_c,Kokkos::LayoutLeft>::value,
                      "C does not have Kokkos::LayoutLeft.");
        
        static_assert(std::is_same<value_type,value_type_c>::value,
                      "A and C do not have the same value type.");
        
        const ordinal_type 
          n = C.extent(0), 
          k = (std::is_same<ArgTrans,Trans::NoTranspose>::value ? A.extent(1) : A.extent(0));
        if (n > 0 && k > 0) {
          if (get_team_rank(member) == 0) {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
            Teuchos::BLAS<ordinal_type,value_type> blas;
            
            blas.HERK(static_cast<const Teuchos::EUplo>(ArgUplo::teuchos),
                      static_cast<const Teuchos::ETransp>(ArgTrans::teuchos),
                      n, k,
                      alpha,
                      A.data(), A.stride_1(),
                      beta,
                      C.data(), C.stride_1());
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
