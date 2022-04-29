#ifndef __KOKKOSBATCHED_SVD_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_SVD_SERIAL_IMPL_HPP__

/// \author Brian Kelley (bmkelle@sandia.gov)

#include "KokkosBatched_SVD_Serial_Internal.hpp"

namespace KokkosBatched {
// Version which computes the full factorization
template <typename AViewType, typename UViewType, typename VViewType,
          typename SViewType, typename WViewType>
KOKKOS_INLINE_FUNCTION int SerialSVD::invoke(SVD_USV_Tag, const AViewType &A,
                                             const UViewType &U,
                                             const SViewType &sigma,
                                             const VViewType &Vt,
                                             const WViewType &work) {
  using value_type = typename AViewType::non_const_value_type;
  return KokkosBatched::SerialSVDInternal::invoke<value_type>(
      A.extent(0), A.extent(1), A.data(), A.stride(0), A.stride(1), U.data(),
      U.stride(0), U.stride(1), Vt.data(), Vt.stride(0), Vt.stride(1),
      sigma.data(), sigma.stride(0), work.data());
}

// Version which computes only singular values
template <typename AViewType, typename SViewType, typename WViewType>
KOKKOS_INLINE_FUNCTION int SerialSVD::invoke(SVD_S_Tag, const AViewType &A,
                                             const SViewType &sigma,
                                             const WViewType &work) {
  using value_type = typename AViewType::non_const_value_type;
  return KokkosBatched::SerialSVDInternal::invoke<value_type>(
      A.extent(0), A.extent(1), A.data(), A.stride(0), A.stride(1), nullptr, 0,
      0, nullptr, 0, 0, sigma.data(), sigma.stride(0), work.data());
}

}  // namespace KokkosBatched

#endif
