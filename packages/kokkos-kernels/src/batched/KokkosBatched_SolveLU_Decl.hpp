#ifndef __KOKKOSBATCHED_SOLVELU_DECL_HPP__
#define __KOKKOSBATCHED_SOLVELU_DECL_HPP__


/// \author Vinh Dang (vqdang@sandia.gov)

#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {
  namespace Experimental {
      
    template<typename ArgAlgo,
             typename TransType>
    struct SerialSolveLU {
      // no piv version
      template<typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const AViewType &A,
             const BViewType &B);
    };       

    template<typename MemberType,
             typename ArgAlgo,
             typename TransType>
    struct TeamSolveLU {
      // no piv version
      template<typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const MemberType &member, 
             const AViewType &A,
             const BViewType &B);
    };       
      
  }
}

#endif
