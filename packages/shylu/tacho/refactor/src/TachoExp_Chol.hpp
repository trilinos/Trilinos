#ifndef __TACHOEXP_CHOL_HPP__
#define __TACHOEXP_CHOL_HPP__

/// \file TachoExp_Chol.hpp
/// \brief Front interface for Cholesky dense factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

namespace Tacho {
  
  namespace Experimental {
    
    ///
    /// Chol:
    /// 
    /// 
    // template<typename SchedulerType,
    //          typename MemoryPoolType,
    //          typename MemberType,
    //          typename ViewTypeA>
    // KOKKOS_INLINE_FUNCTION
    // static int invoke(const SchedulerType &sched,
    //                   const MemoryPoolType &pool,
    //                   const MemberType &member,
    //                   const ViewTypeA &A) {
    //   printf(">> Template Args - Uplo %d, Algo %d\n", 
    //          ArgUplo::tag, ArgAlgo::tag);
    //   TACHO_TEST_FOR_ABORT( true, MSG_INVALID_TEMPLATE_ARGS );
    //   return -1;
    // }

    template<typename ArgUplo, 
             typename ArgAlgo>
    struct Chol;
  }
}

#endif
