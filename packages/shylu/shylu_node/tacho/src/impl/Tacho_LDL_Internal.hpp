#ifndef __TACHO_LDL_INTERNAL_HPP__
#define __TACHO_LDL_INTERNAL_HPP__

/// \file  Tacho_LDL_Internal.hpp
/// \brief LDL team factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Lapack_Team.hpp"

namespace Tacho {

  /// LAPACK LDL
  /// ==========
  template<>
  struct LDL<Uplo::Lower,Algo::Internal> {
    template<typename MemberType,
             typename ViewTypeA,
             typename ViewTypeP,
             typename ViewTypeD,
             typename ViewTypeW>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(MemberType &member,
           const ViewTypeA &A,
           const ViewTypeP &P,
           const ViewTypeD &D,
           const ViewTypeW &W) {
      //typedef typename ViewTypeA::non_const_value_type value_type;
      //typedef typename ViewTypeP::non_const_value_type p_value_type;
        
      static_assert(ViewTypeA::rank == 2,"A is not rank 2 view.");
      static_assert(ViewTypeP::rank == 1,"P is not rank 1 view.");
      static_assert(ViewTypeD::rank == 2,"D is not rank 2 view.");
      static_assert(ViewTypeW::rank == 1,"W is not rank 1 view.");

      /// not done yet
      int r_val(0);
      return r_val;
    }

  };

}

#endif
