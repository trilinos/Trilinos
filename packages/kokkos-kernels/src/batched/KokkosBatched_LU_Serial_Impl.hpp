#ifndef __KOKKOSBATCHED_LU_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_LU_SERIAL_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_LU_Serial_Internal.hpp"


namespace KokkosBatched {
  namespace Experimental {
    ///
    /// Serial Impl
    /// =========

    ///
    /// SerialLU no piv
    ///

#if                                                     \
  defined(__KOKKOSBATCHED_INTEL_MKL__) &&               \
  defined(__KOKKOSBATCHED_INTEL_MKL_BATCHED__) &&       \
  defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
    template<>
    template<typename AViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialLU<Algo::LU::CompactMKL>::
    invoke(const AViewType &A,
           const typename MagnitudeScalarType<typename AViewType::non_const_value_type>::type tiny) {
      typedef typename AViewType::value_type vector_type;
      typedef typename vector_type::value_type value_type;

      const int
        m = A.dimension(0),
        n = A.dimension(1),
        vl = vector_type::vector_length;

      int r_val = 0;
      if (A.stride_0() == 1) {
        LAPACKE_dgetrf_compact(CblasColMajor, 
                               m, n, 
                               (double*)A.data(), A.stride_1(), 
                               (MKL_INT)vl, (MKL_INT)1);
      } else if (A.stride_1() == 1) {
        LAPACKE_dgetrf_compact(CblasRowMajor, 
                               m, n, 
                               (double*)A.data(), A.stride_0(), 
                               (MKL_INT)vl, (MKL_INT)1);
      } else {
        r_val = -1;
      }
      return r_val;
    }
#endif

    template<>
    template<typename AViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialLU<Algo::LU::Unblocked>::
    invoke(const AViewType &A,
           const typename MagnitudeScalarType<typename AViewType::non_const_value_type>::type tiny) {
      return SerialLU_Internal<Algo::LU::Unblocked>::invoke(A.extent(0), A.extent(1),
                                                            A.data(), A.stride_0(), A.stride_1(),
                                                            tiny);
    }
    
    template<>
    template<typename AViewType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialLU<Algo::LU::Blocked>::
    invoke(const AViewType &A,
           const typename MagnitudeScalarType<typename AViewType::non_const_value_type>::type tiny) {
      return SerialLU_Internal<Algo::LU::Blocked>::invoke(A.extent(0), A.extent(1),
                                                          A.data(), A.stride_0(), A.stride_1(),
                                                          tiny);
    }

  }
}

#endif
