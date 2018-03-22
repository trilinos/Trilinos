#ifndef __KOKKOSBATCHED_SCALE_IMPL_HPP__
#define __KOKKOSBATCHED_SCALE_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Scale_Internal.hpp"


namespace KokkosBatched {
  namespace Experimental {
    ///
    /// Serial Impl
    /// ===========
    template<typename ScalarType,
             typename AViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialScale::
    invoke(const ScalarType alpha,
           const AViewType &A) {
      return SerialScaleInternal::
        invoke(A.extent(0), A.extent(1),
               alpha, 
               A.data(), A.stride_0(), A.stride_1());
    }

    ///
    /// Team Impl
    /// =========
    
    template<typename MemberType>
    template<typename ScalarType,
             typename AViewType>
    KOKKOS_INLINE_FUNCTION
    int
    TeamScale<MemberType>::
    invoke(const MemberType &member, 
           const ScalarType alpha,
           const AViewType &A) {
      return TeamScaleInternal::
        invoke(member, 
               A.extent(0), A.extent(1),
               alpha, 
               A.data(), A.stride_0(), A.stride_1());
    }
  
  }
}


#endif
