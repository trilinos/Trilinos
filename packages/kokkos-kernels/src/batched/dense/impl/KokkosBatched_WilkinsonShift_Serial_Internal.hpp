#ifndef __KOKKOSBATCHED_WILKINSON_SHIFT_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_WILKINSON_SHIFT_SERIAL_INTERNAL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {
  ///
  /// Serial Internal Impl
  /// ==================== 
  ///
  /// this impl follows the flame interface of householder transformation
  ///
  struct SerialWilkinsonShiftInternal {
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const ValueType a, const ValueType b,
           const ValueType c, const ValueType d,
           /* */ Kokkos::complex<ValueType> * lambda1,   
           /* */ Kokkos::complex<ValueType> * lambda2,   
           /* */ bool * is_complex) {
      /// compute eigenvalues of 2x2 system [a b;
      ///                                    c d]
      /// when the system has a real complex values,
      /// lambda1 and lambda2 are real eigenvalues
      /// if the system has a complex eigenvalue pair, 
      /// then lambda1 and lambda2 are redefined as follows
      ///   lambda1 := lambda1 + lambda2
      ///   lambda2 := lambda1 * lambda2
      typedef ValueType value_type;
          
      const value_type half(0.5);
      const value_type p = (a+d)*half;
      const value_type q = (b*c-a*d);
      const value_type v = p*p+q;

      if (v < 0) {
        // complex 
        const value_type sqrt_v = Kokkos::Details::ArithTraits<value_type>::sqrt(-v);
        *lambda1 = Kokkos::complex<value_type>(p, sqrt_v);
        *lambda2 = Kokkos::complex<value_type>(p,-sqrt_v);
        *is_complex = true;
      } else {
        // real
        const value_type sqrt_v = Kokkos::Details::ArithTraits<value_type>::sqrt(v);
        *lambda1 = Kokkos::complex<value_type>(p+sqrt_v);
        *lambda2 = Kokkos::complex<value_type>(p-sqrt_v);
        *is_complex = false;
      }
      return 0;
    }
  };

} // end namespace KokkosBatched


#endif
