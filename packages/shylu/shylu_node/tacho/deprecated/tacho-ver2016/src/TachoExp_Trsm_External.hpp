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
        typedef typename ViewTypeA::array_layout array_layout_a;
        typedef typename ViewTypeB::array_layout array_layout_b;
        
        static_assert(ViewTypeA::rank == 2,"A is not rank 2 view.");
        static_assert(ViewTypeB::rank == 2,"B is not rank 2 view.");
        
        static_assert(std::is_same<array_layout_a,Kokkos::LayoutLeft>::value,
                      "A does not have Kokkos::LayoutLeft.");
        static_assert(std::is_same<array_layout_b,Kokkos::LayoutLeft>::value,
                      "B does not have Kokkos::LayoutLeft.");
        
        static_assert(std::is_same<value_type,value_type_b>::value,
                      "A and B do not have the same value type.");
        
        const ordinal_type m = B.extent(0), n = B.extent(1);
        
        if (m > 0 && n > 0) {
          if (get_team_rank(member) == 0) {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
            Teuchos::BLAS<ordinal_type,value_type> blas;
            blas.TRSM(static_cast<const Teuchos::ESide>(ArgSide::teuchos), 
                      static_cast<const Teuchos::EUplo>(ArgUplo::teuchos), 
                      static_cast<const Teuchos::ETransp>(ArgTransA::teuchos), 
                      static_cast<const Teuchos::EDiag>(diagA.teuchos),
                      m, n,
                      alpha,
                      A.data(), A.stride_1(),
                      B.data(), B.stride_1());
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
