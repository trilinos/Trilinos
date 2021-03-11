#ifndef __TACHO_LDL_EXTERNAL_HPP__
#define __TACHO_LDL_EXTERNAL_HPP__

/// \file  Tacho_LDL_External.hpp
/// \brief LAPACK LDL factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Lapack_External.hpp"

namespace Tacho {

  /// LAPACK LDL
  /// ==========
  template<>
  struct LDL<Uplo::Lower,Algo::External> {
    template<typename ViewTypeA,
             typename ViewTypeP,
             typename ViewTypeD,
             typename ViewTypeW>
    inline
    static int
    invoke(const ViewTypeA &A,
           const ViewTypeP &P,
           const ViewTypeD &D,
           const ViewTypeW &W) {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
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
        Lapack<value_type>::sytrf('L',
                                  m,
                                  A.data(), A.stride_1(),
                                  P.data(),
                                  W.data(), W.extent(0),
                                  &r_val);
          
        TACHO_TEST_FOR_EXCEPTION(r_val, std::runtime_error,
                                 "LAPACK (sytrf) returns non-zero error code.");
          
        /// extract diag
        {
          const value_type one(1), zero(0);
          for (ordinal_type i=0;i<m;) {
            const ordinal_type piv = P(i);
            if (piv > 0) {
              /// 1x1 block 
              D(i,0) = A(i,i);
              D(i,1) = zero;
              A(i,i) = one;
              ++i;
            } else if (piv < 0) {
              /// 2x2 symmetric block
              D(i,  0) = A(i,  i  );
              D(i+1,1) = A(i+1,i+1);

              /// this offdiag value should be non-zero
              const value_type offdiag = A(i+1,i);
              D(i,  1) = offdiag;
              D(i+1,0) = offdiag;
              
              A(i,i) = one;
              A(i+1,i+1) = one;
              A(i+1,i) = zero;
            }
          }
        }
      }
#else
      TACHO_TEST_FOR_ABORT( true, ">> This function is only allowed in host space." );
#endif
      return r_val;
    }

  };

}

#endif
