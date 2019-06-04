#ifndef __KOKKOSBATCHED_VECTOR_AVX512D_HPP__
#define __KOKKOSBATCHED_VECTOR_AVX512D_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)
#if defined(__AVX512F__)

#include <immintrin.h>

namespace KokkosBatched {
  namespace Experimental {
    ///
    /// AVX512D double
    ///
    template<>
    class Vector<AVX<double>,8> {
    public:
      using type = Vector<AVX<double>,8>;
      using value_type = double;
      using mag_type = double;

      static const int vector_length = 8;

      union data_type {
        __m512d v;
        double d[8];
      };

      KOKKOS_INLINE_FUNCTION
      static const char* label() { return "AVX512"; }

    private:
      mutable data_type _data;

    public:
      inline Vector() { _data.v = _mm512_setzero_pd(); }
      inline Vector(const value_type val) { _data.v = _mm512_set1_pd(val); }
      inline Vector(const type &b) { _data.v = b._data.v; }
      inline Vector(__m512d const &val) { _data.v = val; }

      inline
      type& operator=(__m512d const &val) {
        _data.v = val;
        return *this;
      }

      inline
      operator __m512d() const {
        return _data.v;
      }

      inline
      type& loadAligned(value_type const *p) {
        _data.v = _mm512_load_pd(p);
        return *this;
      }

      inline
      type& loadUnaligned(value_type const *p) {
        _data.v = _mm512_loadu_pd(p);
        return *this;
      }

      inline
      void storeAligned(value_type *p) const {
        _mm512_store_pd(p, _data.v);
      }

      inline
      void storeUnaligned(value_type *p) const {
        _mm512_storeu_pd(p, _data.v);
      }

      inline
      value_type& operator[](int i) const {
        return _data.d[i];
      }
    };

    
    inline
    static Vector<AVX<double>,8> 
    operator + (const Vector<AVX<double>,8> &a, const Vector<AVX<double>,8> &b) {
      return _mm512_add_pd(a, b);
    }

    
    inline
    static Vector<AVX<double>,8> 
    operator + (const Vector<AVX<double>,8> &a, const double b) {
      return a + Vector<AVX<double>,8>(b);
    }

    
    inline
    static Vector<AVX<double>,8> 
    operator + (const double a, const Vector<AVX<double>,8> &b) {
      return Vector<AVX<double>,8>(a) + b;
    }

    
    inline
    static Vector<AVX<double>,8> & 
    operator += (Vector<AVX<double>,8> &a, const Vector<AVX<double>,8> &b) {
      a = a + b;
      return a;
    }

    
    inline
    static Vector<AVX<double>,8> & 
    operator += (Vector<AVX<double>,8> &a, const double b) {
      a = a + b;
      return a;
    }

    
    inline
    static Vector<AVX<double>,8> 
    operator ++ (Vector<AVX<double>,8> &a, int) {
      Vector<AVX<double>,8> a0 = a;
      a = a + 1.0;
      return a0;
    }

    
    inline
    static Vector<AVX<double>,8> & 
    operator ++ (Vector<AVX<double>,8> &a) {
      a = a + 1.0;
      return a;
    }

    
    inline
    static Vector<AVX<double>,8> 
    operator - (const Vector<AVX<double>,8> &a, const Vector<AVX<double>,8> &b) {
      return _mm512_sub_pd(a, b);
    }

    
    inline
    static Vector<AVX<double>,8> 
    operator - (const Vector<AVX<double>,8> &a, const double b) {
      return a - Vector<AVX<double>,8>(b);
    }

    
    inline
    static Vector<AVX<double>,8> 
    operator - (const double a, const Vector<AVX<double>,8> &b) {
      return Vector<AVX<double>,8>(a) - b;
    }

    
    inline
    static Vector<AVX<double>,8> & 
    operator -= (Vector<AVX<double>,8> &a, const Vector<AVX<double>,8> &b) {
      a = a - b;
      return a;
    }

    
    inline
    static Vector<AVX<double>,8> & 
    operator -= (Vector<AVX<double>,8> &a, const double b) {
      a = a - b;
      return a;
    }

    
    inline
    static Vector<AVX<double>,8> 
    operator -- (Vector<AVX<double>,8> &a, int) {
      Vector<AVX<double>,8> a0 = a;
      a = a - 1.0;
      return a0;
    }

    
    inline
    static Vector<AVX<double>,8> & 
    operator -- (Vector<AVX<double>,8> &a) {
      a = a - 1.0;
      return a;
    }

    
    inline
    static Vector<AVX<double>,8> 
    operator * (const Vector<AVX<double>,8> &a, const Vector<AVX<double>,8> &b) {
      return _mm512_mul_pd(a, b);
    }

    
    inline
    static Vector<AVX<double>,8> 
    operator * (const Vector<AVX<double>,8> &a, const double b) {
      return a * Vector<AVX<double>,8>(b);
    }

    
    inline
    static Vector<AVX<double>,8> 
    operator * (const double a, const Vector<AVX<double>,8> &b) {
      return Vector<AVX<double>,8>(a) * b;
    }

    
    inline
    static Vector<AVX<double>,8> & 
    operator *= (Vector<AVX<double>,8> &a, const Vector<AVX<double>,8> &b) {
      a = a * b;
      return a;
    }

    
    inline
    static Vector<AVX<double>,8> & 
    operator *= (Vector<AVX<double>,8> &a, const double b) {
      a = a * b;
      return a;
    }

    
    inline
    static Vector<AVX<double>,8> 
    operator / (const Vector<AVX<double>,8> &a, const Vector<AVX<double>,8> &b) {
      return _mm512_div_pd(a, b);
    }

    
    inline
    static Vector<AVX<double>,8> 
    operator / (const Vector<AVX<double>,8> &a, const double b) {
      return a / Vector<AVX<double>,8>(b);
    }

    
    inline
    static Vector<AVX<double>,8> 
    operator / (const double a, const Vector<AVX<double>,8> &b) {
      return Vector<AVX<double>,8>(a) / b;
    }

    
    inline
    static Vector<AVX<double>,8> & 
    operator /= (Vector<AVX<double>,8> &a, const Vector<AVX<double>,8> &b) {
      a = a / b;
      return a;
    }

    
    inline
    static Vector<AVX<double>,8> & 
    operator /= (Vector<AVX<double>,8> &a, const double b) {
      a = a / b;
      return a;
    }

    
    inline
    static Vector<AVX<double>,8> 
    operator - (const Vector<AVX<double>,8> &a) {
      return -1*a;
    }

    ///
    /// AVX512D Kokkos::complex<double>
    ///


    template<>
    class Vector<AVX<Kokkos::complex<double> >,4> {
    public:
      using type = Vector<AVX<double>,4> >;
      using value_type = Kokkos::complex<double>;
      using mag_type = double;

      static const int vector_length = 4;

      union data_type {
        __m512d v;
        Kokkos::complex<double> d[4];

        data_type() { v = _mm512_setzero_pd(); }
        data_type(const data_type &b) { v = b.v; }
      };

      KOKKOS_INLINE_FUNCTION
      static const char* label() { return "AVX512"; }

    private:
      mutable data_type _data;

    public:
      inline Vector() { _data.v = _mm512_setzero_pd(); }
      inline Vector(const value_type val) {
        _data.v = _mm512_mask_broadcast_f64x4(_mm512_set1_pd(val.imag()), 0x55, _mm256_set1_pd(val.real()));
      }
      inline Vector(const mag_type val) {
        _data.v = _mm512_mask_broadcast_f64x4(_mm512_setzero_pd(), 0x55, _mm256_set1_pd(val));
      }
      inline Vector(const type &b) { _data.v = b._data.v; }
      inline Vector(__m512d const &val) { _data.v = val; }

      inline
      type& operator=(__m512d const &val) {
        _data.v = val;
        return *this;
      }

      inline
      operator __m512d() const {
        return _data.v;
      }

      inline
      type& loadAligned(value_type const *p) {
        _data.v = _mm512_load_pd((mag_type*)p);
        return *this;
      }

      inline
      type& loadUnaligned(value_type const *p) {
        _data.v = _mm512_loadu_pd((mag_type*)p);
        return *this;
      }

      inline
      void storeAligned(value_type *p) const {
        _mm512_store_pd((mag_type*)p, _data.v);
      }

      inline
      void storeUnaligned(value_type *p) const {
        _mm512_storeu_pd((mag_type*)p, _data.v);
      }

      inline
      value_type& operator[](int i) const {
        return _data.d[i];
      }

    };

    /// complex vector , complex vector

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator + (const Vector<AVX<Kokkos::complex<double> >,4> &a, const Vector<AVX<Kokkos::complex<double> >,4> &b) {
      return _mm512_add_pd(a, b);
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> & 
    operator += (Vector<AVX<Kokkos::complex<double> >,4> &a, const Vector<AVX<Kokkos::complex<double> >,4> &b) {
      a = a + b;
      return a;
    }

    /// complex vector , complex

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator + (const Vector<AVX<Kokkos::complex<double> >,4> &a, const Kokkos::complex<double> b) {
      return a + Vector<AVX<Kokkos::complex<double> >,4>(b);
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator + (const Kokkos::complex<double> a, const Vector<AVX<Kokkos::complex<double> >,4> &b) {
      return Vector<AVX<Kokkos::complex<double> >,4>(a) + b;
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> & 
    operator += (Vector<AVX<Kokkos::complex<double> >,4> &a, const Kokkos::complex<double> b) {
      a = a + b;
      return a;
    }

    /// complex vector , real

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator + (const Vector<AVX<Kokkos::complex<double> >,4> &a, const double b) {
      return a + Vector<AVX<Kokkos::complex<double> >,4>(b);
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator + (const double a, const Vector<AVX<Kokkos::complex<double> >,4> &b) {
      return Vector<AVX<Kokkos::complex<double> >,4>(a) + b;
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> & 
    operator += (Vector<AVX<Kokkos::complex<double> >,4> &a, const double b) {
      a = a + b;
      return a;
    }

    /// unitary operator

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator ++ (Vector<AVX<Kokkos::complex<double> >,4> &a, int) {
      Vector<AVX<Kokkos::complex<double> >,4> a0 = a;
      a = a + 1.0;
      return a0;
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> & 
    operator ++ (Vector<AVX<Kokkos::complex<double> >,4> &a) {
      a = a + 1.0;
      return a;
    }

    /// complex vector , complex vector

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator - (const Vector<AVX<Kokkos::complex<double> >,4> &a, const Vector<AVX<Kokkos::complex<double> >,4> &b) {
      return _mm512_sub_pd(a, b);
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> & 
    operator -= (Vector<AVX<Kokkos::complex<double> >,4> &a, const Vector<AVX<Kokkos::complex<double> >,4> &b) {
      a = a - b;
      return a;
    }

    /// unitary operator

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator -- (Vector<AVX<Kokkos::complex<double> >,4> &a, int) {
      Vector<AVX<Kokkos::complex<double> >,4> a0 = a;
      a = a - 1.0;
      return a0;
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> & 
    operator -- (Vector<AVX<Kokkos::complex<double> >,4> &a) {
      a = a - 1.0;
      return a;
    }

    /// complex vector , complex

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator - (const Vector<AVX<Kokkos::complex<double> >,4> &a, const Kokkos::complex<double> b) {
      return a - Vector<AVX<Kokkos::complex<double> >,4>(b);
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator - (const Kokkos::complex<double> a, const Vector<AVX<Kokkos::complex<double> >,4> &b) {
      return Vector<AVX<Kokkos::complex<double> >,4>(a) - b;
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> & 
    operator -= (Vector<AVX<Kokkos::complex<double> >,4> &a, const Kokkos::complex<double> b) {
      a = a - b;
      return a;
    }

    /// complex vector , real

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator - (const Vector<AVX<Kokkos::complex<double> >,4> &a, const double b) {
      return a - Vector<AVX<Kokkos::complex<double> >,4>(b);
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator - (const double a, const Vector<AVX<Kokkos::complex<double> >,4> &b) {
      return Vector<AVX<Kokkos::complex<double> >,4>(a) - b;
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> & 
    operator -= (Vector<AVX<Kokkos::complex<double> >,4> &a, const double b) {
      a = a - b;
      return a;
    }

    /// complex vector , complex vector

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator * (const Vector<AVX<Kokkos::complex<double> >,4> &a, const Vector<AVX<Kokkos::complex<double> >,4> &b) {
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

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> & 
    operator *= (Vector<AVX<Kokkos::complex<double> >,4> &a, const Vector<AVX<Kokkos::complex<double> >,4> &b) {
      a = a * b;
      return a;
    }

    /// complex vector , complex

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator * (const Vector<AVX<Kokkos::complex<double> >,4> &a, const Kokkos::complex<double> b) {
      return a * Vector<AVX<Kokkos::complex<double> >,4>(b);
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator * (const Kokkos::complex<double> a, const Vector<AVX<Kokkos::complex<double> >,4> &b) {
      return Vector<AVX<Kokkos::complex<double> >,4>(a) * b;
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> & 
    operator *= (Vector<AVX<Kokkos::complex<double> >,4> &a, const Kokkos::complex<double> b) {
      a = a * b;
      return a;
    }

    /// complex vector , real

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator * (const Vector<AVX<Kokkos::complex<double> >,4> &a, const double b) {
      return _mm512_mul_pd(a, _mm512_set1_pd(b));
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator * (const double a, const Vector<AVX<Kokkos::complex<double> >,4> &b) {
      return _mm512_mul_pd(_mm512_set1_pd(a), b);
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> & 
    operator *= (Vector<AVX<Kokkos::complex<double> >,4> &a, const double b) {
      a = a * b;
      return a;
    }

    /// complex vector , complex vector

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator / (const Vector<AVX<Kokkos::complex<double> >,4> &a, const Vector<AVX<Kokkos::complex<double> >,4> &b) {
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

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> & 
    operator /= (Vector<AVX<Kokkos::complex<double> >,4> &a, const Vector<AVX<Kokkos::complex<double> >,4> &b) {
      a = a / b;
      return a;
    }

    /// complex vector , complex

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator / (const Vector<AVX<Kokkos::complex<double> >,4> &a, const Kokkos::complex<double> b) {
      return a / Vector<AVX<Kokkos::complex<double> >,4>(b);
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator / (const Kokkos::complex<double> a, const Vector<AVX<Kokkos::complex<double> >,4> &b) {
      return Vector<AVX<Kokkos::complex<double> >,4>(a) / b;
    }


    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> & 
    operator /= (Vector<AVX<Kokkos::complex<double> >,4> &a, const Kokkos::complex<double> b) {
      a = a / b;
      return a;
    }

    /// complex vector , real

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator / (const Vector<AVX<Kokkos::complex<double> >,4> &a, const double b) {
      return _mm512_div_pd(a, _mm512_set1_pd(b));
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator / (const double a, const Vector<AVX<Kokkos::complex<double> >,4> &b) {
      return Vector<AVX<Kokkos::complex<double> >,4>(a) / b;
    }

  
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> & 
    operator /= (Vector<AVX<Kokkos::complex<double> >,4> &a, const double b) {
      a = a / b;
      return a;
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,4> 
    operator - (const Vector<AVX<Kokkos::complex<double> >,4> &a) {
      return -1.0*a;
    }

  }
}
#endif
#endif
