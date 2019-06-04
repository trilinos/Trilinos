#ifndef __KOKKOSBATCHED_ADD_RADIAL_IMPL_HPP__
#define __KOKKOSBATCHED_ADD_RADIAL_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_AddRadial_Internal.hpp"


namespace KokkosBatched {
  namespace Experimental {
    ///
    /// Serial Impl
    /// ===========
      
    template<typename ScalarType,
             typename AViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialAddRadial::
    invoke(const ScalarType tiny,
           const AViewType &A) {
      return SerialAddRadialInternal::
        invoke((A.extent(0) < A.extent(1) ? A.extent(0) : A.extent(1)),
               tiny, 
               A.data(), (A.stride_0() + A.stride_1()));
    }

    ///
    /// Team Impl
    /// =========
              
    template<typename MemberType>
    template<typename ScalarType,
             typename AViewType>
    KOKKOS_INLINE_FUNCTION
    int
    TeamAddRadial<MemberType>::
    invoke(const MemberType &member, 
           const ScalarType tiny,
           const AViewType &A) {
      return TeamAddRadialInternal::
        invoke(member, 
               (A.extent(0) < A.extent(1) ? A.extent(0) : A.extent(1)),
               tiny, 
               A.data(), (A.stride_0() + A.stride_1()));
    }
    
  } // end namespace Experimental
} //end namespace KokkosBatched


#endif
