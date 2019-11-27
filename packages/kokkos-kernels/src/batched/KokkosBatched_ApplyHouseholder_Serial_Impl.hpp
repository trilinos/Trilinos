#ifndef __KOKKOSBATCHED_APPLY_HOUSEHOLDER_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_APPLY_HOUSEHOLDER_SERIAL_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Householder_Serial_Internal.hpp"

namespace KokkosBatched {

  ///
  /// Serial Impl
  /// ===========

  template<>
  template<typename uViewType,
           typename tauViewType,
           typename AViewType,
           typename wViewType>
  KOKKOS_INLINE_FUNCTION
  int
  SerialApplyHouseholder<Side::Left>::
  invoke(const uViewType &u2,
         const tauViewType &tau,
         const AViewType &A,
         const wViewType &w) {
    return SerialApplyLeftHouseholderInternal::
      invoke(A.extent(0)-1, A.extent(1),
             tau.data(),
             u2.data(), u2.stride(0),
             A.data(), A.stride(1), 
             A.data()+A.stride(0), A.stride(0), A.stride(1),
             w.data());
  }

  template<>
  template<typename uViewType,
           typename tauViewType,
           typename AViewType,
           typename wViewType>
  KOKKOS_INLINE_FUNCTION
  int
  SerialApplyHouseholder<Side::Right>::
  invoke(const uViewType &u2,
         const tauViewType &tau,
         const AViewType &A,
         const wViewType &w) {
    return SerialApplyRightHouseholderInternal::
      invoke(A.extent(0), A.extent(1)-1,
             tau.data(),
             u2.data(), u2.stride(0),
             A.data(), A.stride(0), 
             A.data()+A.stride(1), A.stride(0), A.stride(1),
             w.data());
  }
        
}


#endif
