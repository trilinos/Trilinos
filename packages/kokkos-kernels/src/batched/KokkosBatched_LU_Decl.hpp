#ifndef __KOKKOSBATCHED_LU_DECL_HPP__
#define __KOKKOSBATCHED_LU_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {
  namespace Experimental {
      
    template<typename ArgAlgo>
    struct SerialLU {
      // no piv version
      template<typename AViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const AViewType &A,
             const typename MagnitudeScalarType<typename AViewType::non_const_value_type>::type tiny = 0);
    };       

    template<typename MemberType,
             typename ArgAlgo>
    struct TeamLU {
      // no piv version
      template<typename AViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const MemberType &member, 
             const AViewType &A,
             const typename MagnitudeScalarType<typename AViewType::non_const_value_type>::type tiny = 0);
    };       
      
  }
}

#endif
