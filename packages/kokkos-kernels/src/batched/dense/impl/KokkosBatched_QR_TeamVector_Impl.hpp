#ifndef __KOKKOSBATCHED_QR_TEAMVECTOR_IMPL_HPP__
#define __KOKKOSBATCHED_QR_TEAMVECTOR_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_QR_TeamVector_Internal.hpp"

namespace KokkosBatched {

  ///
  /// TeamVector Impl
  /// ===============

  template<typename MemberType>
  struct TeamVectorQR<MemberType,Algo::QR::Unblocked> {
    template<typename AViewType,
             typename tViewType,
             typename wViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member,
           const AViewType &A,
           const tViewType &t,
           const wViewType &w) {
      return TeamVectorQR_Internal::
        invoke(member,
               A.extent(0), A.extent(1), 
               A.data(), A.stride_0(), A.stride_1(),
               t.data(), t.stride_0(), 
               w.data());
    }
  };
        
}



#endif
