#ifndef __TACHOEXP_HERK_HPP__
#define __TACHOEXP_HERK_HPP__

/// \file TachoExp_Herk.hpp
/// \brief Front interface for Herk operators
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

namespace Tacho {

  namespace Experimental {
    
    template<typename ArgUplo, typename ArgTrans, typename ArgAlgo>
    struct Herk {
      template<typename PolicyType,
               typename MemberType,
               typename ScalarType,
               typename ViewTypeA,
               typename ViewTypeC>
      KOKKOS_INLINE_FUNCTION
      static int invoke(const PolicyType &policy,
                        const MemberType &member,
                        const ScalarType alpha,
                        const ViewTypeA &A,
                        const ScalarType beta,
                        const ViewTypeC &C) { 
        printf(">> Template Args - Uplo %d, Trans %d, Algo %d\n", 
               ArgUplo::tag, ArgTrans::tag, ArgAlgo::tag);  
        TACHO_TEST_FOR_ABORT( true, MSG_INVALID_TEMPLATE_ARGS );
        return -1;
      }
    };
  }
}

#endif
