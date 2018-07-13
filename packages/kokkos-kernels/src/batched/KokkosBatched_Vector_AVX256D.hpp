#ifndef __KOKKOSBATCHED_VECTOR_AVX256D_HPP__
#define __KOKKOSBATCHED_VECTOR_AVX256D_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#if defined(__AVX__) || defined(__AVX2__)

#include <immintrin.h>

namespace KokkosBatched {
  namespace Experimental {

    ///
    /// AVX256D double
    ///

    template<>
    class Vector<AVX<double>,4> {
    public:
      using type = Vector<AVX<double>,4>;
      using value_type = double;
      using mag_type = double;

      static const int vector_length = 4;

      union data_type {
        __m256d v;
        double d[4];
      };

      KOKKOS_INLINE_FUNCTION
      static const char* label() { return "AVX256"; }

    private:
      mutable data_type _data;

    public:
      inline Vector() { _data.v = _mm256_setzero_pd(); }
      inline Vector(const value_type val) { _data.v = _mm256_set1_pd(val); }
      inline Vector(const type &b) { _data.v = b._data.v; }
      inline Vector(__m256d const &val) { _data.v = val; }

      inline
      type& operator=(__m256d const &val) {
        _data.v = val;
        return *this;
      }

      inline
      operator __m256d() const {
        return _data.v;
      }

      inline
      type& loadAligned(value_type const *p) {
        _data.v = _mm256_load_pd(p);
        return *this;
      }

      inline
      type& loadUnaligned(value_type const *p) {
        _data.v = _mm256_loadu_pd(p);
        return *this;
      }

      inline
      void storeAligned(value_type *p) const {
        _mm256_store_pd(p, _data.v);
      }

      inline
      void storeUnaligned(value_type *p) const {
        _mm256_storeu_pd(p, _data.v);
      }

      inline
      value_type& operator[](int i) const {
        return _data.d[i];
      }
    };

    
    inline
    static Vector<AVX<double>,4> 
    operator + (Vector<AVX<double>,4> const & a, Vector<AVX<double>,4> const & b) {
      return _mm256_add_pd(a, b);
    }

    
    inline
    static Vector<AVX<double>,4> 
    operator + (Vector<AVX<double>,4> const & a, const double b) {
      return a + Vector<AVX<double>,4>(b);
    }

    
    inline
    static Vector<AVX<double>,4> 
    operator + (const double a, Vector<AVX<double>,4> const & b) {
      return Vector<AVX<double>,4>(a) + b;
    }

    
    inline
    static Vector<AVX<double>,4> & 
    operator += (Vector<AVX<double>,4> & a, Vector<AVX<double>,4> const & b) {
      a = a + b;
      return a;
    }

    
    inline
    static Vector<AVX<double>,4> & 
    operator += (Vector<AVX<double>,4> & a, const double b) {
      a = a + b;
      return a;
    }

    
    inline
    static Vector<AVX<double>,4> 
    operator ++ (Vector<AVX<double>,4> & a, int) {
      Vector<AVX<double>,4> a0 = a;
      a = a + 1.0;
      return a0;
    }

    
    inline
    static Vector<AVX<double>,4> & 
    operator ++ (Vector<AVX<double>,4> & a) {
      a = a + 1.0;
      return a;
    }

    
    inline
    static Vector<AVX<double>,4> 
    operator - (Vector<AVX<double>,4> const & a, Vector<AVX<double>,4> const & b) {
      return _mm256_sub_pd(a, b);
    }

    
    inline
    static Vector<AVX<double>,4> 
    operator - (Vector<AVX<double>,4> const & a, const double b) {
      return a - Vector<AVX<double>,4>(b);
    }

    
    inline
    static Vector<AVX<double>,4> 
    operator - (const double a, Vector<AVX<double>,4> const & b) {
      return Vector<AVX<double>,4>(a) - b;
    }

    
    inline
    static Vector<AVX<double>,4> & 
    operator -= (Vector<AVX<double>,4> & a, Vector<AVX<double>,4> const & b) {
      a = a - b;
      return a;
    }

    
    inline
    static Vector<AVX<double>,4> & 
    operator -= (Vector<AVX<double>,4> & a, const double b) {
      a = a - b;
      return a;
    }

    
    inline
    static Vector<AVX<double>,4> 
    operator -- (Vector<AVX<double>,4> & a, int) {
      Vector<AVX<double>,4> a0 = a;
      a = a - 1.0;
      return a0;
    }

    
    inline
    static Vector<AVX<double>,4> & 
    operator -- (Vector<AVX<double>,4> & a) {
      a = a - 1.0;
      return a;
    }

    
    inline
    static Vector<AVX<double>,4> 
    operator * (Vector<AVX<double>,4> const & a, Vector<AVX<double>,4> const & b) {
      return _mm256_mul_pd(a, b);
    }

    
    inline
    static Vector<AVX<double>,4> 
    operator * (Vector<AVX<double>,4> const & a, const double b) {
      return a * Vector<AVX<double>,4>(b);
    }

    
    inline
    static Vector<AVX<double>,4> operator * (const double a, Vector<AVX<double>,4> const & b) {
      return Vector<AVX<double>,4>(a) * b;
    }

    
    inline
    static Vector<AVX<double>,4> & 
    operator *= (Vector<AVX<double>,4> & a, Vector<AVX<double>,4> const & b) {
      a = a * b;
      return a;
    }

    
    inline
    static Vector<AVX<double>,4> & 
    operator *= (Vector<AVX<double>,4> & a, const double b) {
      a = a * b;
      return a;
    }

    
    inline
    static Vector<AVX<double>,4> 
    operator / (Vector<AVX<double>,4> const & a, Vector<AVX<double>,4> const & b) {
      return _mm256_div_pd(a, b);
    }

    
    inline
    static Vector<AVX<double>,4> 
    operator / (Vector<AVX<double>,4> const & a, const double b) {
      return a / Vector<AVX<double>,4>(b);
    }

    
    inline
    static Vector<AVX<double>,4> 
    operator / (const double a, Vector<AVX<double>,4> const & b) {
      return Vector<AVX<double>,4>(a) / b;
    }

    
    inline
    static Vector<AVX<double>,4> & 
    operator /= (Vector<AVX<double>,4> & a, Vector<AVX<double>,4> const & b) {
      a = a / b;
      return a;
    }

    
    inline
    static Vector<AVX<double>,4> & 
    operator /= (Vector<AVX<double>,4> & a, const double b) {
      a = a / b;
      return a;
    }

    
    inline
    static Vector<AVX<double>,4> 
    operator - (Vector<AVX<double>,4> const & a) {
      return -1.0*a;
    }

    ///
    /// AVX256D Kokkos::complex<double>
    //

    template<>
    class Vector<AVX<Kokkos::complex<double> >,2> {
    public:
      using type = Vector<AVX<Kokkos::complex<double> >,2>;
      using value_type = Kokkos::complex<double>;
      using mag_type = double;

      static const int vector_length = 2;

      union data_type {
        __m256d v;
        Kokkos::complex<double> d[2];

        data_type() { v = _mm256_setzero_pd(); }
        data_type(const data_type &b) { v = b.v; }
      };

      KOKKOS_INLINE_FUNCTION
      static const char* label() { return "AVX256"; }

    private:
      mutable data_type _data;

    public:
      inline Vector() { _data.v = _mm256_setzero_pd(); }
      inline Vector(const value_type val) { _data.v = _mm256_broadcast_pd((__m128d const *)&val);}
      inline Vector(const mag_type val) { const value_type a(val); _data.v = _mm256_broadcast_pd((__m128d const *)&a); }
      inline Vector(const type &b) { _data.v = b._data.v; }
      inline Vector(__m256d const &val) { _data.v = val; }

      inline
      type& operator=(__m256d const &val) {
        _data.v = val;
        return *this;
      }

      inline
      operator __m256d() const {
        return _data.v;
      }

      inline
      type& loadAligned(value_type const *p) {
        _data.v = _mm256_load_pd((mag_type*)p);
        return *this;
      }

      inline
      type& loadUnaligned(value_type const *p) {
        _data.v = _mm256_loadu_pd((mag_type*)p);
        return *this;
      }

      inline
      void storeAligned(value_type *p) const {
        _mm256_store_pd((mag_type*)p, _data.v);
      }

      inline
      void storeUnaligned(value_type *p) const {
        _mm256_storeu_pd((mag_type*)p, _data.v);
      }

      inline
      value_type& operator[](int i) const {
        return _data.d[i];
      }
    };

    /// complex vector , complex vector

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator + (Vector<AVX<Kokkos::complex<double> >,2> const & a, Vector<AVX<Kokkos::complex<double> >,2> const & b) {
      return _mm256_add_pd(a, b);
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> & 
    operator += (Vector<AVX<Kokkos::complex<double> >,2> & a, Vector<AVX<Kokkos::complex<double> >,2> const & b) {
      a = a + b;
      return a;
    }

    /// complex vector , complex

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator + (Vector<AVX<Kokkos::complex<double> >,2> const & a, const Kokkos::complex<double> b) {
      return a + Vector<AVX<Kokkos::complex<double> >,2>(b);
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator + (const Kokkos::complex<double> a, Vector<AVX<Kokkos::complex<double> >,2> const & b) {
      return Vector<AVX<Kokkos::complex<double> >,2>(a) + b;
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> & 
    operator += (Vector<AVX<Kokkos::complex<double> >,2> & a, const Kokkos::complex<double> b) {
      a = a + b;
      return a;
    }

    /// complex vector , real

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator + (Vector<AVX<Kokkos::complex<double> >,2> const & a, const double b) {
      return a + Vector<AVX<Kokkos::complex<double> >,2>(b);
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator + (const double a, Vector<AVX<Kokkos::complex<double> >,2> const & b) {
      return Vector<AVX<Kokkos::complex<double> >,2>(a) + b;
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> & 
    operator += (Vector<AVX<Kokkos::complex<double> >,2> & a, const double b) {
      a = a + b;
      return a;
    }

    /// unitary operator

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator ++ (Vector<AVX<Kokkos::complex<double> >,2> & a, int) {
      Vector<AVX<Kokkos::complex<double> >,2> a0 = a;
      a = a + 1.0;
      return a0;
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> & 
    operator ++ (Vector<AVX<Kokkos::complex<double> >,2> & a) {
      a = a + 1.0;
      return a;
    }

    /// complex vector , complex vector

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator - (Vector<AVX<Kokkos::complex<double> >,2> const & a, Vector<AVX<Kokkos::complex<double> >,2> const & b) {
      return _mm256_sub_pd(a, b);
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> & 
    operator -= (Vector<AVX<Kokkos::complex<double> >,2> & a, Vector<AVX<Kokkos::complex<double> >,2> const & b) {
      a = a - b;
      return a;
    }

    /// complex vector , complex

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator - (Vector<AVX<Kokkos::complex<double> >,2> const & a, const Kokkos::complex<double> b) {
      return a - Vector<AVX<Kokkos::complex<double> >,2>(b);
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator - (const Kokkos::complex<double> a, Vector<AVX<Kokkos::complex<double> >,2> const & b) {
      return Vector<AVX<Kokkos::complex<double> >,2>(a) - b;
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> & 
    operator -= (Vector<AVX<Kokkos::complex<double> >,2> & a, const Kokkos::complex<double> b) {
      a = a - b;
      return a;
    }

    /// complex vector , real

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator - (Vector<AVX<Kokkos::complex<double> >,2> const & a, const double b) {
      return a - Vector<AVX<Kokkos::complex<double> >,2>(b);
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator - (const double a, Vector<AVX<Kokkos::complex<double> >,2> const & b) {
      return Vector<AVX<Kokkos::complex<double> >,2>(a) - b;
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> & 
    operator -= (Vector<AVX<Kokkos::complex<double> >,2> & a, const double b) {
      a = a - b;
      return a;
    }

    /// unitary operator

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator -- (Vector<AVX<Kokkos::complex<double> >,2> & a, int) {
      Vector<AVX<Kokkos::complex<double> >,2> a0 = a;
      a = a - 1.0;
      return a0;
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> & 
    operator -- (Vector<AVX<Kokkos::complex<double> >,2> & a) {
      a = a - 1.0;
      return a;
    }

    /// complex vector , complex vector

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator * (Vector<AVX<Kokkos::complex<double> >,2> const & a, Vector<AVX<Kokkos::complex<double> >,2> const & b) {
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
    
    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> & 
    operator *= (Vector<AVX<Kokkos::complex<double> >,2> & a, Vector<AVX<Kokkos::complex<double> >,2> const & b) {
      a = a * b;
      return a;
    }

    /// complex vector , complex

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator * (Vector<AVX<Kokkos::complex<double> >,2> const & a, const Kokkos::complex<double> b) {
      return a * Vector<AVX<Kokkos::complex<double> >,2>(b);
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator * (const Kokkos::complex<double> a, Vector<AVX<Kokkos::complex<double> >,2> const & b) {
      return Vector<AVX<Kokkos::complex<double> >,2>(a) * b;
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> & 
    operator *= (Vector<AVX<Kokkos::complex<double> >,2> & a, const Kokkos::complex<double> b) {
      a = a * b;
      return a;
    }

    /// complex vector , real

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator * (Vector<AVX<Kokkos::complex<double> >,2> const & a, const double b) {
      return _mm256_mul_pd(a, _mm256_set1_pd(b));
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator * (const double a, Vector<AVX<Kokkos::complex<double> >,2> const & b) {
      return _mm256_mul_pd(_mm256_set1_pd(a), b);
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> & 
    operator *= (Vector<AVX<Kokkos::complex<double> >,2> & a, const double b) {
      a = a * b;
      return a;
    }

    /// complex vector , complex vector

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator / (Vector<AVX<Kokkos::complex<double> >,2> const & a, Vector<AVX<Kokkos::complex<double> >,2> const & b) {
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

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> & 
    operator /= (Vector<AVX<Kokkos::complex<double> >,2> & a, Vector<AVX<Kokkos::complex<double> >,2> const & b) {
      a = a / b;
      return a;
    }

    /// complex vector , complex

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator / (Vector<AVX<Kokkos::complex<double> >,2> const & a, const Kokkos::complex<double> b) {
      return a / Vector<AVX<Kokkos::complex<double> >,2>(b);
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator / (const Kokkos::complex<double> a, Vector<AVX<Kokkos::complex<double> >,2> const & b) {
      return Vector<AVX<Kokkos::complex<double> >,2>(a) / b;
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> & 
    operator /= (Vector<AVX<Kokkos::complex<double> >,2> & a, const Kokkos::complex<double> b) {
      a = a / b;
      return a;
    }

    /// complex vector , real

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator / (Vector<AVX<Kokkos::complex<double> >,2> const & a, const double b) {
      return _mm256_div_pd(a, _mm256_set1_pd(b));
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator / (const double a, Vector<AVX<Kokkos::complex<double> >,2> const & b) {
      return Vector<AVX<Kokkos::complex<double> >,2>(a) / b;
    }

    
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> & 
    operator /= (Vector<AVX<Kokkos::complex<double> >,2> & a, const double b) {
      a = a / b;
      return a;
    }

   
    inline
    static Vector<AVX<Kokkos::complex<double> >,2> 
    operator - (Vector<AVX<Kokkos::complex<double> >,2> const & a) {
      return -1.0*a;
    }

  }
}

#endif
#endif
