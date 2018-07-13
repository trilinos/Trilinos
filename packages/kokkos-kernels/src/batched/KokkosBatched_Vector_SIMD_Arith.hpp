#ifndef __KOKKOSBATCHED_VECTOR_SIMD_ARITH_HPP__
#define __KOKKOSBATCHED_VECTOR_SIMD_ARITH_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Kokkos_Complex.hpp"

namespace KokkosBatched {
  namespace Experimental {

    //#define KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(A) typename std::enable_if<std::is_same<Kokkos::Impl::ActiveExecutionMemorySpace,Kokkos::HostSpace>::value,A >::type
#define KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(T,l) Vector< SIMD< T >, l >
#define KOKKOSKERNELS_SIMD_ARITH_RETURN_REFERENCE_TYPE(T,l) Vector< SIMD < T >, l> &

    /// simd, simd
#if defined(__KOKKOSBATCHED_ENABLE_AVX__)
#if defined(__AVX512F__)
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(double,8)
    operator + (const Vector<SIMD<double>,8> &a, const Vector<SIMD<double>,8> &b) {
      return _mm512_add_pd(a, b);
    }

#if !defined(KOKKOS_COMPILER_GNU)
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<double>,4)
    operator + (const Vector<SIMD<Kokkos::complex<double> >,4> &a, const Vector<SIMD<Kokkos::complex<double> >,4> &b) {
      return _mm512_add_pd(a, b);
    }
#endif

#endif
#if defined(__AVX__) || defined(__AVX2__)
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(double,4)
    operator + (const Vector<SIMD<double>,4> & a, const Vector<SIMD<double>,4> & b) {
      return _mm256_add_pd(a, b);
    }

#if !defined(KOKKOS_COMPILER_GNU)
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<double>,2)
    operator + (const Vector<SIMD<Kokkos::complex<double> >,2> & a, const Vector<SIMD<Kokkos::complex<double> >,2> & b) {
      return _mm256_add_pd(a, b);
    }
#endif

#endif
#endif

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION 
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(T,l)
    operator + (const Vector<SIMD<T>,l> &a,  const Vector<SIMD<T>,l> &b) {
      Vector<SIMD<T>,l> r_val;
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
      for (int i=0;i<l;++i)
        r_val[i] = a[i] + b[i];
      return r_val;
    }
        
    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_REFERENCE_TYPE(T,l)
    operator += (Vector<SIMD<T>,l> &a, const Vector<SIMD<T>,l> &b) {
      a = a + b;
      return a;
    }
    
    /// simd, real
    
    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(T,l)
    operator + (const Vector<SIMD<T>,l> &a, const T b) {
      return a + Vector<SIMD<T>,l>(b);
    }
    
    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(T,l)
    operator + (const T a, const Vector<SIMD<T>,l> &b) {
      return Vector<SIMD<T>,l>(a) + b;
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_REFERENCE_TYPE(T,l)
    operator += (Vector<SIMD<T>,l> &a, const T b) {
      a = a + b;
      return a;
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(T,l)
    operator ++ (Vector<SIMD<T>,l> &a, int) {
      Vector<SIMD<T>,l> a0 = a;
      a = a + typename Kokkos::Details::ArithTraits<T>::mag_type(1);
      return a0;
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_REFERENCE_TYPE(T,l)
    operator ++ (Vector<SIMD<T>,l> &a) {
      a = a + typename Kokkos::Details::ArithTraits<T>::mag_type(1);
      return a;
    }

    /// simd complex, real

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<T>,l)
    operator + (const Vector<SIMD<Kokkos::complex<T> >,l> &a, const T b) {
      return a + Vector<SIMD<Kokkos::complex<T> >,l>(b);
    }
    
    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<T>,l)
    operator + (const T a, const Vector<SIMD<Kokkos::complex<T> >,l> &b) {
      return Vector<SIMD<Kokkos::complex<T> >,l>(a) + b;
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_REFERENCE_TYPE(Kokkos::complex<T>,l)
    operator += (Vector<SIMD<Kokkos::complex<T> >,l> &a, const T b) {
      a = a + b;
      return a;
    }

    /// simd complex, complex 

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<T>,l)
    operator + (const Vector<SIMD<Kokkos::complex<T> >,l> &a, const Kokkos::complex<T> b) {
      return a + Vector<SIMD<Kokkos::complex<T> >,l>(b);
    }
    
    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<T>,l)
    operator + (const Kokkos::complex<T> a, const Vector<SIMD<Kokkos::complex<T> >,l> &b) {
      return Vector<SIMD<Kokkos::complex<T> >,l>(a) + b;
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_REFERENCE_TYPE(Kokkos::complex<T>,l)
    operator += (Vector<SIMD<Kokkos::complex<T> >,l> &a, const Kokkos::complex<T> b) {
      a = a + b;
      return a;
    }

    /// ---------------------------------------------------------------------------------------------------

    /// simd, simd

#if defined(__KOKKOSBATCHED_ENABLE_AVX__)
#if defined(__AVX512F__)
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(double,8)
    operator - (const Vector<SIMD<double>,8> &a, const Vector<SIMD<double>,8> &b) {
      return _mm512_sub_pd(a, b);
    }

#if !defined(KOKKOS_COMPILER_GNU)
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<double>,4)
    operator - (const Vector<SIMD<Kokkos::complex<double> >,4> &a, const Vector<SIMD<Kokkos::complex<double> >,4> &b) {
      return _mm512_sub_pd(a, b);
    }
#endif

#endif
#if defined(__AVX__) || defined(__AVX2__)
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(double,4)
    operator - (const Vector<SIMD<double>,4> & a, const Vector<SIMD<double>,4> & b) {
      return _mm256_sub_pd(a, b);
    }

#if !defined(KOKKOS_COMPILER_GNU)
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<double>,2)
    operator - (const Vector<SIMD<Kokkos::complex<double> >,2> & a, const Vector<SIMD<Kokkos::complex<double> >,2> & b) {
      return _mm256_sub_pd(a, b);
    }
#endif

#endif
#endif

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(T,l)
    operator - (const Vector<SIMD<T>,l> &a, const Vector<SIMD<T>,l> &b) {
      Vector<SIMD<T>,l> r_val;
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
      for (int i=0;i<l;++i)
        r_val[i] = a[i] - b[i];
      return r_val;
    }


    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(T,l)
    operator - (const Vector<SIMD<T>,l> &a) {
      Vector<SIMD<T>,l> r_val;
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
      for (int i=0;i<l;++i)
        r_val[i] = -a[i];
      return r_val;
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_REFERENCE_TYPE(T,l)
    operator -= (Vector<SIMD<T>,l> &a, const Vector<SIMD<T>,l> &b) {
      a = a - b;
      return a;
    }

    /// simd, real

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(T,l)
    operator - (const Vector<SIMD<T>,l> &a, const T b) {
      return a - Vector<SIMD<T>,l>(b);
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(T,l)
    operator - (const T a, const Vector<SIMD<T>,l> &b) {
      return Vector<SIMD<T>,l>(a) - b;
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_REFERENCE_TYPE(T,l)
    operator -= (Vector<SIMD<T>,l> &a, const T b) {
      a = a - b;
      return a;
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(T,l)
    operator -- (Vector<SIMD<T>,l> &a, int) {
      Vector<SIMD<T>,l> a0 = a;
      a = a - typename Kokkos::Details::ArithTraits<T>::mag_type(1);
      return a0;
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_REFERENCE_TYPE(T,l)
    operator -- (Vector<SIMD<T>,l> &a) {
      a = a - typename Kokkos::Details::ArithTraits<T>::mag_type(1);
      return a;
    }

    /// simd complex, real
    
    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<T>,l)
    operator - (const Vector<SIMD<Kokkos::complex<T> >,l> &a, const T b) {
      return a - Vector<SIMD<Kokkos::complex<T> >,l>(b);
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<T>,l)
    operator - (const T a, const Vector<SIMD<Kokkos::complex<T> >,l> &b) {
      return Vector<SIMD<Kokkos::complex<T> >,l>(a) - b;
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_REFERENCE_TYPE(Kokkos::complex<T>,l)
    operator -= (Vector<SIMD<Kokkos::complex<T> >,l> &a, const T b) {
      a = a - b;
      return a;
    }

    /// simd complex, complex
    
    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<T>,l)
    operator - (const Vector<SIMD<Kokkos::complex<T> >,l> &a, const Kokkos::complex<T> b) {
      return a - Vector<SIMD<Kokkos::complex<T> >,l>(b);
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<T>,l)
    operator - (const Kokkos::complex<T> a, const Vector<SIMD<Kokkos::complex<T> >,l> &b) {
      return Vector<SIMD<Kokkos::complex<T> >,l>(a) - b;
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_REFERENCE_TYPE(Kokkos::complex<T>,l)
    operator -= (Vector<SIMD<Kokkos::complex<T> >,l> &a, const Kokkos::complex<T> b) {
      a = a - b;
      return a;
    }

    /// ---------------------------------------------------------------------------------------------------    

    /// simd, simd

#if defined(__KOKKOSBATCHED_ENABLE_AVX__)
#if defined(__AVX512F__)
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(double,8)
    operator * (const Vector<SIMD<double>,8> &a, const Vector<SIMD<double>,8> &b) {
      return _mm512_mul_pd(a, b);
    }

#if !defined(KOKKOS_COMPILER_GNU)
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<double>,4)
    operator * (const Vector<SIMD<Kokkos::complex<double> >,4> &a, const Vector<SIMD<Kokkos::complex<double> >,4> &b) {
      const __m512d
        as = _mm512_permute_pd(a, 0x55),
        br = _mm512_permute_pd(b, 0x00),
        bi = _mm512_permute_pd(b, 0xff);
      
#if defined(__FMA__)
      // latency 7, throughput 0.5
      return _mm512_fmaddsub_pd(a, br, _mm512_mul_pd(as, bi));
#else
      return _mm512_add_pd(_mm512_mul_pd(a, br),
                           _mm512_castsi512_pd(_mm512_xor_si512(_mm512_castpd_si512(_mm512_mul_pd(as, bi)),
                                                                _mm512_castpd_si512(_mm512_mask_broadcast_f64x4(_mm512_setzero_pd(), 0x55,
                                                                                                                _mm256_set1_pd(-0.0))))));
      // const __mm512d cc = _mm512_mul_pd(as, bi);
      // return _mm512_mask_sub_pd(_mm512_mask_add_pd(_mm512_mul_pd(a, br), 0x55, cc), 0xaa, cc);
#endif
    }
#endif

#endif
#if defined(__AVX__) || defined(__AVX2__)
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(double,4)
    operator * (const Vector<SIMD<double>,4> & a, const Vector<SIMD<double>,4> & b) {
      return _mm256_mul_pd(a, b);
    }

#if !defined(KOKKOS_COMPILER_GNU)
    KOKKOS_FORCEINLINE_FUNCTION
    static Vector<SIMD<Kokkos::complex<double> >,2>
    operator * (const Vector<SIMD<Kokkos::complex<double> >,2> & a, const Vector<SIMD<Kokkos::complex<double> >,2> & b) {
      const __m256d
        as = _mm256_permute_pd(a, 0x5),
        br = _mm256_permute_pd(b, 0x0),
        bi = _mm256_permute_pd(b, 0xf);
      
#if defined(__FMA__)
      return _mm256_fmaddsub_pd(a, br, _mm256_mul_pd(as, bi));
#else
      return _mm256_add_pd(_mm256_mul_pd(a, br),
                           _mm256_xor_pd(_mm256_mul_pd(as, bi),
                                         _mm256_set_pd( 0.0, -0.0, 0.0, -0.0)));
#endif
    }
#endif

#endif
#endif

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(T,l)
    operator * (const Vector<SIMD<T>,l> &a, const Vector<SIMD<T>,l> &b) {
      Vector<SIMD<T>,l> r_val;
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
      for (int i=0;i<l;++i)
        r_val[i] = a[i] * b[i];
      return r_val;
    }
        
    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_REFERENCE_TYPE(T,l)
    operator *= (Vector<SIMD<T>,l> &a, const Vector<SIMD<T>,l> &b) {
      a = a * b;
      return a;
    }


    /// simd, real

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(T,l)
    operator * (const Vector<SIMD<T>,l> &a, const T b) {
      return a * Vector<SIMD<T>,l>(b);
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(T,l)
    operator * (const T a, const Vector<SIMD<T>,l> &b) {
      return Vector<SIMD<T>,l>(a) * b;
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_REFERENCE_TYPE(T,l)
    operator *= (Vector<SIMD<T>,l> &a, const T b) {
      a = a * b;
      return a;
    }

    /// simd complex, real

#if defined(__KOKKOSBATCHED_ENABLE_AVX__)
#if defined(__AVX512F__)

#if !defined(KOKKOS_COMPILER_GNU)
    KOKKOS_FORCEINLINE_FUNCTION
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<double>,4)
    operator * (const Vector<SIMD<Kokkos::complex<double> >,4> &a, const double b) {
      return _mm512_mul_pd(a, _mm512_set1_pd(b));
    }
#endif

#endif
#if defined(__AVX__) || defined(__AVX2__)

#if !defined(KOKKOS_COMPILER_GNU)
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<double>,2)
    operator * (const Vector<SIMD<Kokkos::complex<double> >,2> & a, const double b) {
      return _mm256_mul_pd(a, _mm256_set1_pd(b));
    }
#endif

#endif
#endif

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<T>,l)
    operator * (const Vector<SIMD<Kokkos::complex<T> >,l> &a, const T b) {
      return a * Vector<SIMD<Kokkos::complex<T> >,l>(b);
    }

#if defined(__KOKKOSBATCHED_ENABLE_AVX__)
#if defined(__AVX512F__)

#if !defined(KOKKOS_COMPILER_GNU)
    KOKKOS_FORCEINLINE_FUNCTION
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<double>,4)
    operator * (const double a, const Vector<SIMD<Kokkos::complex<double> >,4> &b) {
      return _mm512_mul_pd(_mm512_set1_pd(a), b);
    }
#endif

#endif
#if defined(__AVX__) || defined(__AVX2__)

#if !defined(KOKKOS_COMPILER_GNU)
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<double>,2)
    operator * (const double a, const Vector<SIMD<Kokkos::complex<double> >,2> & b) {
      return _mm256_mul_pd(_mm256_set1_pd(a), b);
    }
#endif

#endif
#endif

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<T>,l)
    operator * (const T a, const Vector<SIMD<Kokkos::complex<T> >,l> &b) {
      return Vector<SIMD<Kokkos::complex<T> >,l>(a) * b;
    }


    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_REFERENCE_TYPE(Kokkos::complex<T>,l)
    operator *= (Vector<SIMD<Kokkos::complex<T> >,l> &a, const T b) {
      a = a * b;
      return a;
    }

    /// simd complex, complex

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<T>,l)
    operator * (const Vector<SIMD<Kokkos::complex<T> >,l> &a, const Kokkos::complex<T> b) {
      return a * Vector<SIMD<Kokkos::complex<T> >,l>(b);
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<T>,l)
    operator * (const Kokkos::complex<T> a, const Vector<SIMD<Kokkos::complex<T> >,l> &b) {
      return Vector<SIMD<Kokkos::complex<T> >,l>(a) * b;
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_REFERENCE_TYPE(Kokkos::complex<T>,l)
    operator *= (Vector<SIMD<Kokkos::complex<T> >,l> &a, const Kokkos::complex<T> b) {
      a = a * b;
      return a;
    }

    /// ---------------------------------------------------------------------------------------------------    

    /// simd, simd

#if defined(__KOKKOSBATCHED_ENABLE_AVX__)
#if defined(__AVX512F__)
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(double,8)
    operator / (const Vector<SIMD<double>,8> &a, const Vector<SIMD<double>,8> &b) {
      return _mm512_div_pd(a, b);
    }

#if !defined(KOKKOS_COMPILER_GNU)
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<double>,4)
    operator / (const Vector<SIMD<Kokkos::complex<double> >,4> &a, const Vector<SIMD<Kokkos::complex<double> >,4> &b) {
      const __m512d
        as = _mm512_permute_pd(a, 0x55),
        cb = _mm512_castsi512_pd(_mm512_xor_si512(_mm512_castpd_si512(b),
                                                  _mm512_castpd_si512(_mm512_mask_broadcast_f64x4(_mm512_setzero_pd(), 0xAA,
                                                                                                  _mm256_set1_pd(-0.0))))),
        br = _mm512_permute_pd(cb, 0x00),
        bi = _mm512_permute_pd(cb, 0xff);

#if defined(__FMA__)
      return _mm512_div_pd(_mm512_fmaddsub_pd(a,  br, _mm512_mul_pd(as, bi)),
                           _mm512_fmadd_pd   (br, br, _mm512_mul_pd(bi, bi)));
#else
      return _mm512_div_pd(_mm512_add_pd(_mm512_mul_pd(a, br),
                                         _mm512_castsi512_pd(_mm512_xor_si512(_mm512_castpd_si512(_mm512_mul_pd(as, bi)),
                                                                              _mm512_castpd_si512(_mm512_mask_broadcast_f64x4(_mm512_setzero_pd(), 0xAA,
                                                                                                                              _mm256_set1_pd(-0.0)))))),
                           _mm512_add_pd(_mm512_mul_pd(br, br), _mm512_mul_pd(bi, bi)));
      // const __mm512d cc = _mm512_mul_pd(as, bi);
      // return _mm512_div_pd(_mm512_mask_sub_pd(_mm512_mask_add_pd(_mm512_mul_pd(a, br), 0x55, cc), 0xaa, cc),
      //                      _mm512_add_pd(_mm512_mul_pd(br, br), _mm512_mul_pd(bi, bi)));
#endif
    }
#endif

#endif

#if defined(__AVX__) || defined(__AVX2__)
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(double,4)
    operator / (const Vector<SIMD<double>,4> & a, const Vector<SIMD<double>,4> & b) {
      return _mm256_div_pd(a, b);
    }

#if !defined(KOKKOS_COMPILER_GNU)
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<double>,2)
    operator / (Vector<SIMD<Kokkos::complex<double> >,2> const & a, Vector<SIMD<Kokkos::complex<double> >,2> const & b) {
      const __m256d
        as = _mm256_permute_pd(a, 0x5),
        cb = _mm256_xor_pd(b, _mm256_set_pd(-0.0, 0.0, -0.0, 0.0)),
        br = _mm256_permute_pd(cb, 0x0),
        bi = _mm256_permute_pd(cb, 0xf);

#if defined(__FMA__)
      return _mm256_div_pd(_mm256_fmaddsub_pd(a, br, _mm256_mul_pd(as, bi)),
                           _mm256_add_pd(_mm256_mul_pd(br, br), _mm256_mul_pd(bi, bi)));
#else
      return _mm256_div_pd(_mm256_add_pd(_mm256_mul_pd(a, br),
                                         _mm256_xor_pd(_mm256_mul_pd(as, bi),
                                                       _mm256_set_pd( 0.0, -0.0, 0.0, -0.0))),
                           _mm256_add_pd(_mm256_mul_pd(br, br), _mm256_mul_pd(bi, bi)));
#endif
    }
#endif

#endif
#endif

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(T,l)
    operator / (const Vector<SIMD<T>,l> &a, const Vector<SIMD<T>,l> &b) {
      Vector<SIMD<T>,l> r_val;
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
      for (int i=0;i<l;++i)
        r_val[i] = a[i] / b[i];
      return r_val;
    }


    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_REFERENCE_TYPE(T,l)
    operator /= (Vector<SIMD<T>,l> &a, const Vector<SIMD<T>,l> &b) {
      a = a / b;
      return a;
    }

    /// simd, real
#if defined(__KOKKOSBATCHED_ENABLE_AVX__)
#if defined(__AVX512F__)

#if !defined(KOKKOS_COMPILER_GNU)
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<double>,4)
    operator / (const Vector<SIMD<Kokkos::complex<double> >,4> &a, const double b) {
      return _mm512_div_pd(a, _mm512_set1_pd(b));
    }
#endif

#endif
#endif

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(T,l)
    operator / (const Vector<SIMD<T>,l> &a, const T b) {
      return a / Vector<SIMD<T>,l>(b);
    }


    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(T,l)
    operator / (const T a, const Vector<SIMD<T>,l> &b) {
      return Vector<SIMD<T>,l>(a) / b;
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_REFERENCE_TYPE(T,l)
    operator /= (Vector<SIMD<T>,l> &a, const T b) {
      a = a / b;
      return a;
    }

    /// simd complex, real

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<T>,l)
    operator / (const Vector<SIMD<Kokkos::complex<T> >,l> &a, const T b) {
      return a / Vector<SIMD<Kokkos::complex<T> >,l>(b);
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<T>,l)
    operator / (const T a, const Vector<SIMD<Kokkos::complex<T> >,l> &b) {
      return Vector<SIMD<Kokkos::complex<T> >,l>(a) / b;
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_REFERENCE_TYPE(Kokkos::complex<T>,l)
    operator /= (Vector<SIMD<Kokkos::complex<T> >,l> &a, const T b) {
      a = a / b;
      return a;
    }

    /// simd complex, complex

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static 
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<T>,l)
    operator / (const Vector<SIMD<Kokkos::complex<T> >,l> &a, const Kokkos::complex<T> b) {
      return a / Vector<SIMD<Kokkos::complex<T> >,l>(b);
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE(Kokkos::complex<T>,l)
    operator / (const Kokkos::complex<T> a, const Vector<SIMD<Kokkos::complex<T> >,l> &b) {
      return Vector<SIMD<Kokkos::complex<T> >,l>(a) / b;
    }

    template<typename T, int l>
    KOKKOS_FORCEINLINE_FUNCTION
    static
    KOKKOSKERNELS_SIMD_ARITH_RETURN_REFERENCE_TYPE(Kokkos::complex<T>,l)
    operator /= (Vector<SIMD<Kokkos::complex<T> >,l> &a, const Kokkos::complex<T> b) {
      a = a / b;
      return a;
    }
#undef KOKKOSKERNELS_SIMD_ARITH_RETURN_TYPE
#undef KOKKOSKERNELS_SIMD_ARITH_RETURN_REFERENCE_TYPE
  }
}

#endif
