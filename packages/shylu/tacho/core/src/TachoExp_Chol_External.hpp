#ifndef __TACHOEXP_CHOL_EXTERNAL_HPP__
#define __TACHOEXP_CHOL_EXTERNAL_HPP__

/// \file  TachoExp_Chol_Upper_ExternalLapack.hpp
/// \brief LAPACK upper Cholesky factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Lapack_External.hpp"

namespace Tacho {

  namespace Experimental {

    /// LAPACK Chol
    /// ===========
    template<typename ArgUplo>
    struct Chol<ArgUplo,Algo::External> {
      template<typename SchedType,
               typename MemberType,
               typename ViewTypeA>
      inline
      static int
      invoke(const SchedType &sched,
             const MemberType &member,
             const ViewTypeA &A) {
        typedef typename ViewTypeA::non_const_value_type value_type;
        
        static_assert(ViewTypeA::rank == 2,"A is not rank 2 view.");

        int r_val = 0;      
        
        const ordinal_type m = A.dimension_0();
        if (m > 0) {
          if (get_team_rank(member) == 0) {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
            Lapack<value_type>::potrf(ArgUplo::param,
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
