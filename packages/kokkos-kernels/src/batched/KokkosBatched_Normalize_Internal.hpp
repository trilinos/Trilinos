#ifndef __KOKKOSBATCHED_NORMALIZE_INTERNAL_HPP__
#define __KOKKOSBATCHED_NORMALIZE_INTERNAL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

  ///
  /// Serial Internal Impl
  /// ==================== 
  struct SerialNormalizeInternal {
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const int m,
           /* */ ValueType *__restrict__ v, const int vs) {
      typedef ValueType value_type;
      typedef Kokkos::Details::ArithTraits<value_type> ats;
      typedef typename ats::mag_type mag_type;

      mag_type norm(0);
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int i=0;i<m;++i) {
        const auto v_at_i = v[i*vs];
        norm += ats::real(v_at_i*ats::conj(v_at_i));
      }
      norm = Kokkos::Details::ArithTraits<mag_type>::sqrt(norm);
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int i=0;i<m;++i) 
        v[i*vs] /= norm;

      return 0;
    }

    template<typename RealType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const int m,
           /* */ RealType *__restrict__ vr, const int vrs,
           /* */ RealType *__restrict__ vi, const int vis) {
      typedef RealType real_type;
      typedef Kokkos::Details::ArithTraits<real_type> ats;
      typedef typename ats::mag_type mag_type;

      mag_type norm(0);
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int i=0;i<m;++i) {
        const auto vr_at_i = vr[i*vrs];
        const auto vi_at_i = vi[i*vis];
        norm += vr_at_i*vr_at_i+vi_at_i*vi_at_i;
      }
      norm = Kokkos::Details::ArithTraits<mag_type>::sqrt(norm);
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int i=0;i<m;++i) { 
        vr[i*vrs] /= norm;
        vi[i*vis] /= norm;
      }

      return 0;
    }
  };

} // end namespace KokkosBatched


#endif
