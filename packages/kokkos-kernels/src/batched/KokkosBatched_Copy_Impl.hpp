#ifndef __KOKKOSBATCHED_COPY_IMPL_HPP__
#define __KOKKOSBATCHED_COPY_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Copy_Internal.hpp"


namespace KokkosBatched {
  namespace Experimental {
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
      return CopyInternal::
        invoke(A.dimension_0(), 
               A.dimension_1(), 
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
      return CopyInternal::
        invoke(A.dimension_1(), 
               A.dimension_0(), 
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
        return CopyInternal::
            invoke(member,
                   A.dimension_0(), 
                   A.dimension_1(), 
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
        return CopyInternal::
          invoke(member,
                 A.dimension_1(), 
                 A.dimension_0(), 
                 A.data(), A.stride_1(), A.stride_0(),
                 B.data(), B.stride_0(), B.stride_1());
      }
    };
      
  } // end namespace Experimental
} //end namespace KokkosBatched


#endif
