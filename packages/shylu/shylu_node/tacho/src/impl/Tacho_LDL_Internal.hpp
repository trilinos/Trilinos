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
      typedef typename ViewTypeA::non_const_value_type value_type;
      typedef typename ViewTypeP::non_const_value_type p_value_type;
        
      static_assert(ViewTypeA::rank == 2,"A is not rank 2 view.");
      static_assert(ViewTypeP::rank == 1,"P is not rank 1 view.");
      static_assert(ViewTypeD::rank == 2,"D is not rank 2 view.");
      static_assert(ViewTypeW::rank == 1,"W is not rank 1 view.");

      int r_val = 0;
      const ordinal_type m = A.extent(0);
      if (m > 0) {
        /// factorize LDL
        // LapackTeam<value_type>::sytrf(member,
        //                               Uplo::Lower::param,
        //                               m,
        //                               A.data(), A.stride_1(),
        //                               P.data(),
        //                               W.data(), W.extent(0),
        //                               &r_val);
        // member.team_barrier();

        /// extract diag
        {
          const value_type one(1), zero(0);
          Kokkos::parallel_for
            (Kokkos::TeamVectorRange(member, m),
             [&,P,A,D](const ordinal_type &i) {
              const ordinal_type 
                piv = P(i),
                prev_piv = i > 0 ? P(i-1) : 1;
              if (piv > 0) {
                /// 1x1 block 
                D(i,0) = A(i,i);
                A(i,i) = one;
              } else if (piv < 0 && prev_piv < 0) {
                /// do nothing; 
              } else { /// piv < 0 but prev_piv > 0
                /// 2x2 symmetric block
                D(i,  0) = A(i,  i  );
                D(i+1,1) = A(i+1,i+1);
                
                const value_type offdiag= A(i+1,i);
                D(i,  1) = offdiag;
                D(i+1,0) = offdiag;
                
                A(i,i) = one;
                A(i+1,i+1) = one;
                A(i+1,i) = zero;
              }
            });
        }
      }
      return r_val;
    }

  };

}

#endif
