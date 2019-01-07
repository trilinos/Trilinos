#ifndef __KOKKOSBATCHED_INVERSELU_DECL_HPP__
#define __KOKKOSBATCHED_INVERSELU_DECL_HPP__


/// \author Vinh Dang (vqdang@sandia.gov)

#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {
  namespace Experimental {
      
    template<typename ArgAlgo>
    struct SerialInverseLU {
      // no piv version
      template<typename AViewType,
               typename WViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const AViewType &A,
             const WViewType &W);
    };       

    template<typename MemberType,
             typename ArgAlgo>
    struct TeamInverseLU {
      // no piv version
      template<typename AViewType,
               typename WViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const MemberType &member, 
             const AViewType &A,
             const WViewType &W);
    };       
      
  }
}

#endif
