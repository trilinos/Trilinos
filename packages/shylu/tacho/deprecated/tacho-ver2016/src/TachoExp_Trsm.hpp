#ifndef __TACHOEXP_TRSM_HPP__
#define __TACHOEXP_TRSM_HPP__

/// \file TachoExp_Trsm.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

namespace Tacho {
  
  namespace Experimental {
    template<typename ArgSide,
             typename ArgUplo, 
             typename ArgTrans, 
             typename ArgAlgo>
    struct Trsm {
      template<typename PolicyType,
               typename MemberType,
               typename DiagType,
               typename ScalarType,
               typename ViewTypeA,
               typename ViewTypeB>
      KOKKOS_INLINE_FUNCTION
      static int invoke(const PolicyType &policy,
                        const MemberType &member,
                        const DiagType diagA,
                        const ScalarType alpha,
                        const ViewTypeA &A,
                        const ViewTypeB &B) {
        printf(">> Template Args - Side %d, Uplo %d, Trans %d, Algo %d\n", 
               ArgSide::tag, ArgUplo::tag, ArgTrans::tag, ArgAlgo::tag);
        TACHO_TEST_FOR_ABORT( true, MSG_INVALID_TEMPLATE_ARGS );
        return -1;
      }
    };
  }  
}

#endif
