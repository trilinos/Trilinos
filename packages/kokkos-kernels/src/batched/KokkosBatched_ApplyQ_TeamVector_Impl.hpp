#ifndef __KOKKOSBATCHED_APPLY_Q_TEAMVECTOR_IMPL_HPP__
#define __KOKKOSBATCHED_APPLY_Q_TEAMVECTOR_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_ApplyQ_TeamVector_Internal.hpp"

namespace KokkosBatched {

  ///
  /// TeamVector Impl
  /// ===============

  template<typename MemberType>
  struct TeamVectorApplyQ<MemberType,Side::Left,Trans::NoTranspose,Algo::ApplyQ::Unblocked> {
    template<typename AViewType,
             typename tViewType,
             typename BViewType,
             typename wViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const AViewType &A,
           const tViewType &t,
           const BViewType &B,
           const wViewType &w) {
      return TeamVectorApplyQ_LeftForwardInternal::
        invoke(member,
               B.extent(0), B.extent(1), A.extent(1), 
               A.data(), A.stride_0(), A.stride_1(),
               t.data(), t.stride_0(), 
               B.data(), B.stride_0(), B.stride_1(),
               w.data());
    }
  };

  template<typename MemberType>
  struct TeamVectorApplyQ<MemberType,Side::Left,Trans::Transpose,Algo::ApplyQ::Unblocked> {
    template<typename AViewType,
             typename tViewType,
             typename BViewType,
             typename wViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const AViewType &A,
           const tViewType &t,
           const BViewType &B,
           const wViewType &w) {
      return TeamVectorApplyQ_LeftBackwardInternal::
        invoke(member,
               B.extent(0), B.extent(1), A.extent(1), 
               A.data(), A.stride_0(), A.stride_1(),
               t.data(), t.stride_0(), 
               B.data(), B.stride_0(), B.stride_1(),
               w.data());
    }
  };

  template<typename MemberType>
  struct TeamVectorApplyQ<MemberType,Side::Right,Trans::NoTranspose,Algo::ApplyQ::Unblocked> {
    template<typename AViewType,
             typename tViewType,
             typename BViewType,
             typename wViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const AViewType &A,
           const tViewType &t,
           const BViewType &B,
           const wViewType &w) {
      return TeamVectorApplyQ_RightForwardInternal::
        invoke(member,
               B.extent(0), B.extent(1), A.extent(1), 
               A.data(), A.stride_0(), A.stride_1(),
               t.data(), t.stride_0(), 
               B.data(), B.stride_0(), B.stride_1(),
               w.data());
    }
  };
  
}



#endif
