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
             typename ViewTypeW>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(MemberType &member,
           const ViewTypeA &A,
           const ViewTypeP &P,
           const ViewTypeW &W) {
      typedef typename ViewTypeA::non_const_value_type value_type;
      //typedef typename ViewTypeP::non_const_value_type p_value_type;
        
      static_assert(ViewTypeA::rank == 2,"A is not rank 2 view.");
      static_assert(ViewTypeP::rank == 1,"P is not rank 1 view.");
      static_assert(ViewTypeW::rank == 1,"W is not rank 1 view.");

      TACHO_TEST_FOR_EXCEPTION(P.extent(0) < 4*A.extent(0), std::runtime_error,
                               "P should be 4*A.extent(0) .");

      int r_val(0);
      const ordinal_type m = A.extent(0);
      if (m > 0) {
        /// factorize LDL
        LapackTeam<value_type>::sytrf(member, 
                                      Uplo::Lower::param,
                                      m,
                                      A.data(), A.stride_1(),
                                      P.data()+m, /// fpiv is input
                                      W.data(),
                                      &r_val);
      }
      return r_val;
    }

    template<typename MemberType,
             typename ViewTypeA,
             typename ViewTypeP,
             typename ViewTypeD>
    inline
    static int
    modify(MemberType &member,
           const ViewTypeA &A,
           const ViewTypeP &P,
           const ViewTypeD &D) {
      typedef typename ViewTypeA::non_const_value_type value_type;
        
      static_assert(ViewTypeA::rank == 2,"A is not rank 2 view.");
      static_assert(ViewTypeP::rank == 1,"P is not rank 1 view.");
      static_assert(ViewTypeD::rank == 2,"D is not rank 2 view.");

      TACHO_TEST_FOR_ABORT(D.extent(0) < A.extent(0), 
                           "D extent(0) is smaller than A extent(0).");
      TACHO_TEST_FOR_ABORT(D.extent(1) != 2, 
                           "D is supposed to store 2x2 blocks .");
      TACHO_TEST_FOR_ABORT(P.extent(0) < 4*A.extent(0), 
                           "P should be 4*A.extent(0) .");

      int r_val = 0;      
      const ordinal_type m = A.extent(0);
      if (m > 0) {
        ordinal_type 
          *__restrict__ ipiv = P.data(),
          *__restrict__ fpiv = ipiv + m, 
          *__restrict__ perm = fpiv + m, 
          *__restrict__ peri = perm + m;
        constexpr value_type one(1);
        Kokkos::parallel_for(Kokkos::TeamVectorRange(member,m),[&](const int &i) {
            D(i,0) = A(i,i);
            A(i,i) = one;
            ipiv[i] = i+1;
            fpiv[i] = 0;
            perm[i] = i;
            peri[i] = i;
          });
      }
      return r_val;
    }
    
  };

}

#endif
