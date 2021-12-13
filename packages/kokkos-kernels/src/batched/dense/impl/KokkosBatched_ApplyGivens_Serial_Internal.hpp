#ifndef __KOKKOSBATCHED_APPLY_GIVENS_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_APPLY_GIVENS_SERIAL_INTERNAL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

  ///
  /// Serial Internal Impl
  /// ==================== 
  ///
  /// this impl follows the flame interface of householder transformation
  ///
  struct SerialApplyLeftGivensInternal {
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const Kokkos::pair<ValueType,ValueType> G,
           const int n, 
           /* */ ValueType * a1t, const int a1ts,
           /* */ ValueType * a2t, const int a2ts) {
      typedef ValueType value_type;
      if (n == 0) return 0; // quick return
      if (G.first == value_type(1) && G.second == value_type(0)) return 0;
      /// G = [ gamma -sigma;
      ///       sigma  gamma ];
      /// A := G A
      /// where gamma is G.first and sigma is G.second 
        
      const value_type gamma = G.first;
      const value_type sigma = G.second;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int j=0;j<n;++j) {
        const value_type alpha1 = a1t[j*a1ts];
        const value_type alpha2 = a2t[j*a2ts];
        a1t[j*a1ts] = gamma*alpha1 - sigma*alpha2;
        a2t[j*a1ts] = sigma*alpha1 + gamma*alpha2;
      }
      return 0;
    }
  };

  struct SerialApplyRightGivensInternal {
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const Kokkos::pair<ValueType,ValueType> G,
           const int m, 
           /* */ ValueType * a1, const int a1s,
           /* */ ValueType * a2, const int a2s) {
      typedef ValueType value_type;
      if (m == 0) return 0; // quick return
      if (G.first == value_type(1) && G.second == value_type(0)) return 0;
      /// G = [ gamma -sigma;
      ///       sigma  gamma ];
      /// A := A G'
      /// where gamma is G.first and sigma is G.second 
          
      const value_type gamma = G.first;
      const value_type sigma = G.second;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif        
      for (int i=0;i<m;++i) {
        const value_type alpha1 = a1[i*a1s];
        const value_type alpha2 = a2[i*a2s];
        a1[i*a1s] = gamma*alpha1 - sigma*alpha2;
        a2[i*a1s] = sigma*alpha1 + gamma*alpha2;
      }
      return 0;
    }
  };

  struct SerialApplyLeftRightGivensInternal {
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const Kokkos::pair<ValueType,ValueType> &G12,
           const int &m, const int &n,
           /* */ ValueType *__restrict__ a1t, 
           /* */ ValueType *__restrict__ a2t, 
           /* */ ValueType *__restrict__ a1, 
           /* */ ValueType *__restrict__ a2, 
           const int &as0, const int &as1) {
      typedef ValueType value_type;
      if (G12.first == value_type(1) && G12.second == value_type(0)) return 0;
      if (m == 0 && n == 0) return 0; // quick return

      const value_type gamma12 = G12.first;
      const value_type sigma12 = G12.second;

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int j=0;j<n;++j) {
        const value_type alpha1 = a1t[j*as1];
        const value_type alpha2 = a2t[j*as1];
        a1t[j*as1] = ( gamma12*alpha1 - sigma12*alpha2 );
        a2t[j*as1] = ( sigma12*alpha1 + gamma12*alpha2 );
      }

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int i=0;i<m;++i) {
        const value_type alpha1 = a1[i*as0];
        const value_type alpha2 = a2[i*as0];
        a1[i*as0] = ( gamma12*alpha1 - sigma12*alpha2 );
        a2[i*as0] = ( sigma12*alpha1 + gamma12*alpha2 );
      }
      return 0;
    }

    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const Kokkos::pair<ValueType,ValueType> &G12,
           const Kokkos::pair<ValueType,ValueType> &G13,
           const int &m, const int &n, 
           /* */ ValueType *__restrict__ a1t, 
           /* */ ValueType *__restrict__ a2t, 
           /* */ ValueType *__restrict__ a3t, 
           /* */ ValueType *__restrict__ a1, 
           /* */ ValueType *__restrict__ a2, 
           /* */ ValueType *__restrict__ a3, 
           const int &as0, const int &as1) {
      typedef ValueType value_type;
      if (m == 0 && n == 0) return 0; // quick return

      const value_type gamma12 = G12.first;
      const value_type sigma12 = G12.second;
      const value_type gamma13 = G13.first;
      const value_type sigma13 = G13.second;

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int j=0;j<n;++j) {
        const value_type alpha2 = a2t[j*as1];
        const value_type alpha3 = a3t[j*as1];
        {
          const value_type alpha1 = a1t[j*as1];
          a1t[j*as1] = ( gamma12*alpha1 - sigma12*alpha2 );
          a2t[j*as1] = ( sigma12*alpha1 + gamma12*alpha2 ); 
        }
        {
          const value_type alpha1 = a1t[j*as1];
          a1t[j*as1] = ( gamma13*alpha1 - sigma13*alpha3 );
          a3t[j*as1] = ( sigma13*alpha1 + gamma13*alpha3 );
        }
      }

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int i=0;i<m;++i) {
        const value_type alpha2 = a2[i*as0];
        const value_type alpha3 = a3[i*as0];
        {
          const value_type alpha1 = a1[i*as0];
          a1[i*as0] = ( gamma12*alpha1 - sigma12*alpha2 );
          a2[i*as0] = ( sigma12*alpha1 + gamma12*alpha2 );
        }
        {
          const value_type alpha1 = a1[i*as0];
          a1[i*as0] = ( gamma13*alpha1 - sigma13*alpha3 );
          a3[i*as0] = ( sigma13*alpha1 + gamma13*alpha3 );
        }
      }
      return 0;
    }
  };
    
} // end namespace KokkosBatched


#endif
