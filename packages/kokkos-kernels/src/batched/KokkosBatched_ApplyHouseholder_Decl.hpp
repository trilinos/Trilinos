#ifndef __KOKKOSBATCHED_APPLY_HOUSEHOLDER_DECL_HPP__
#define __KOKKOSBATCHED_APPLY_HOUSEHOLDER_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

  ///
  /// Serial Householder 
  ///

  // level 1 operation
  template<typename ArgSide>
  struct SerialApplyHouseholder {
    template<typename uViewType,
             typename tauViewType,
             typename AViewType,
             typename wViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const uViewType &u2,
           const tauViewType &tau,
           const AViewType
           const wViewType &w);
  };

  // level 1 operation
  template<typename MemberType, 
           typename ArgSide>
  struct TeamVectorApplyHouseholder {
    template<typename uViewType,
             typename tauViewType,
             typename AViewType,
             typename wViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const uViewType &u2,
           const tauViewType &tau,
           const AViewType
           const wViewType &w);
  };

}

#include "KokkosBatched_ApplyHouseholder_Serial_Impl.hpp"
#include "KokkosBatched_ApplyHouseholder_TeamVector_Impl.hpp"

#endif
