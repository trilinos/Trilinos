#ifndef __KOKKOSBATCHED_QR_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_QR_SERIAL_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_QR_Serial_Internal.hpp"

namespace KokkosBatched {

  ///
  /// Serial Impl
  /// ===========

  template<>
  template<typename AViewType,
           typename tViewType,
           typename wViewType>
  KOKKOS_INLINE_FUNCTION
  int
  SerialQR<Algo::QR::Unblocked>::
  invoke(const AViewType &A,
         const tViewType &t,
         const wViewType &w) {
    return SerialQR_Internal::
      invoke(A.extent(0), A.extent(1), 
             A.data(), A.stride_0(), A.stride_1(),
             t.data(), t.stride_0(), 
             w.data());
  }
        
}



#endif
