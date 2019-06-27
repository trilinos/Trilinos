#ifndef __KOKKOSBATCHED_COPY_IMPL_HPP__
#define __KOKKOSBATCHED_COPY_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Copy_Internal.hpp"

namespace KokkosBatched {

  ///
  /// Serial Impl
  /// ===========
      
  template<>
  template<typename AViewType,
           typename BViewType>
  KOKKOS_INLINE_FUNCTION
  int
  SerialCopy<Trans::NoTranspose>::
  invoke(const AViewType &A,
         /* */ BViewType &B) {
    return SerialCopyInternal::
      invoke(A.extent(0), 
             A.extent(1), 
             A.data(), A.stride_0(), A.stride_1(),
             B.data(), B.stride_0(), B.stride_1());
  }
    
  template<>
  template<typename AViewType,
           typename BViewType>
  KOKKOS_INLINE_FUNCTION
  int
  SerialCopy<Trans::Transpose>::
  invoke(const AViewType &A,
         /* */ BViewType &B) {
    return SerialCopyInternal::
      invoke(A.extent(1), 
             A.extent(0), 
             A.data(), A.stride_1(), A.stride_0(),
             B.data(), B.stride_0(), B.stride_1());
  }
    
    
  ///
  /// Team Impl
  /// =========
    
  template<typename MemberType>
  struct TeamCopy<MemberType,Trans::NoTranspose> {
    template<typename AViewType,
             typename BViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const AViewType &A,
           /* */ BViewType &B) {
      return TeamCopyInternal::
        invoke(member,
               A.extent(0), 
               A.extent(1), 
               A.data(), A.stride_0(), A.stride_1(),
               B.data(), B.stride_0(), B.stride_1());
    }
  };
    
  template<typename MemberType>
  struct TeamCopy<MemberType,Trans::Transpose> {
    template<typename AViewType,
             typename BViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const AViewType &A,
           /* */ BViewType &B) {
      return TeamCopyInternal::
        invoke(member,
               A.extent(1), 
               A.extent(0), 
               A.data(), A.stride_1(), A.stride_0(),
               B.data(), B.stride_0(), B.stride_1());
    }
  };
      
} //end namespace KokkosBatched


#endif
