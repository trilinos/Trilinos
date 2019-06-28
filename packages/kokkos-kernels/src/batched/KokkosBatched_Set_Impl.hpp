#ifndef __KOKKOSBATCHED_SET_IMPL_HPP__
#define __KOKKOSBATCHED_SET_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Set_Internal.hpp"

namespace KokkosBatched {

    ///
    /// Serial Impl
    /// ===========
      
    template<typename ScalarType,
             typename AViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialSet::
    invoke(const ScalarType alpha,
           const AViewType &A) {
      return SerialSetInternal::
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
    TeamSet<MemberType>::
    invoke(const MemberType &member, 
           const ScalarType alpha,
           const AViewType &A) {
      return TeamSetInternal::
        invoke(member, 
               A.extent(0), A.extent(1),
               alpha, 
               A.data(), A.stride_0(), A.stride_1());
    }
} //end namespace KokkosBatched


#endif
