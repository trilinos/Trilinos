#ifndef __TACHOEXP_CHOL_EXTERNAL_HPP__
#define __TACHOEXP_CHOL_EXTERNAL_HPP__

/// \file  TachoExp_Chol_Upper_ExternalLapack.hpp
/// \brief LAPACK upper Cholesky factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Teuchos_LAPACK.hpp"

namespace Tacho {

  namespace Experimental {

    /// LAPACK Chol
    /// ===========
    template<typename ArgUplo>
    struct Chol<ArgUplo,Algo::External> {
      template<typename PolicyType,
               typename MemberType,
               typename ViewTypeA>
      inline
      static int
      invoke(const PolicyType &policy,
             const MemberType &member,
             const ViewTypeA &A) {
        typedef typename ViewTypeA::non_const_value_type value_type;
        typedef typename ViewTypeA::array_layout array_layout;
        
        static_assert(ViewTypeA::rank == 2,"A is not rank 2 view.");
        static_assert(std::is_same<array_layout,Kokkos::LayoutLeft>::value, 
                      "A does not have Kokkos::LayoutLeft.");
        int r_val = 0;      
        
        const ordinal_type m = A.extent(0);
        if (m > 0) {
          if (get_team_rank(member) == 0) {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
            Teuchos::LAPACK<ordinal_type,value_type> lapack;

            lapack.POTRF(ArgUplo::param,
                         m, 
                         A.data(), A.stride_1(),
                         &r_val);

            TACHO_TEST_FOR_EXCEPTION(r_val, std::runtime_error, 
                                     "LAPACK (potrf) returns non-zero error code.");
#else
            TACHO_TEST_FOR_ABORT( true, ">> This function is only allowed in host space." );
#endif
          }
        }
        return r_val;
      }
    };
  }
}

#endif
