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

    template<typename SpT>
    class Vector<VectorTag<AVX<double,SpT>,8> > {
    public:
      using type = Vector<VectorTag<AVX<double,SpT>,8> >;
      using value_type = double;
      using real_type = double;

      enum : int { vector_length = 8 };

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

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > operator + (Vector<VectorTag<AVX<double,SpT>,8> > const & a, Vector<VectorTag<AVX<double,SpT>,8> > const & b) {
      return _mm512_add_pd(a, b);
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > operator + (Vector<VectorTag<AVX<double,SpT>,8> > const & a, const double b) {
      return a + Vector<VectorTag<AVX<double,SpT>,8> >(b);
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > operator + (const double a, Vector<VectorTag<AVX<double,SpT>,8> > const & b) {
      return Vector<VectorTag<AVX<double,SpT>,8> >(a) + b;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > & operator += (Vector<VectorTag<AVX<double,SpT>,8> > & a, Vector<VectorTag<AVX<double,SpT>,8> > const & b) {
      a = a + b;
      return a;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > & operator += (Vector<VectorTag<AVX<double,SpT>,8> > & a, const double b) {
      a = a + b;
      return a;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > operator ++ (Vector<VectorTag<AVX<double,SpT>,8> > & a, int) {
      Vector<VectorTag<AVX<double,SpT>,8> > a0 = a;
      a = a + 1.0;
      return a0;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > & operator ++ (Vector<VectorTag<AVX<double,SpT>,8> > & a) {
      a = a + 1.0;
      return a;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > operator - (Vector<VectorTag<AVX<double,SpT>,8> > const & a, Vector<VectorTag<AVX<double,SpT>,8> > const & b) {
      return _mm512_sub_pd(a, b);
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > operator - (Vector<VectorTag<AVX<double,SpT>,8> > const & a, const double b) {
      return a - Vector<VectorTag<AVX<double,SpT>,8> >(b);
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > operator - (const double a, Vector<VectorTag<AVX<double,SpT>,8> > const & b) {
      return Vector<VectorTag<AVX<double,SpT>,8> >(a) - b;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > & operator -= (Vector<VectorTag<AVX<double,SpT>,8> > & a, Vector<VectorTag<AVX<double,SpT>,8> > const & b) {
      a = a - b;
      return a;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > & operator -= (Vector<VectorTag<AVX<double,SpT>,8> > & a, const double b) {
      a = a - b;
      return a;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > operator -- (Vector<VectorTag<AVX<double,SpT>,8> > & a, int) {
      Vector<VectorTag<AVX<double,SpT>,8> > a0 = a;
      a = a - 1.0;
      return a0;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > & operator -- (Vector<VectorTag<AVX<double,SpT>,8> > & a) {
      a = a - 1.0;
      return a;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > operator * (Vector<VectorTag<AVX<double,SpT>,8> > const & a, Vector<VectorTag<AVX<double,SpT>,8> > const & b) {
      return _mm512_mul_pd(a, b);
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > operator * (Vector<VectorTag<AVX<double,SpT>,8> > const & a, const double b) {
      return a * Vector<VectorTag<AVX<double,SpT>,8> >(b);
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > operator * (const double a, Vector<VectorTag<AVX<double,SpT>,8> > const & b) {
      return Vector<VectorTag<AVX<double,SpT>,8> >(a) * b;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > & operator *= (Vector<VectorTag<AVX<double,SpT>,8> > & a, Vector<VectorTag<AVX<double,SpT>,8> > const & b) {
      a = a * b;
      return a;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > & operator *= (Vector<VectorTag<AVX<double,SpT>,8> > & a, const double b) {
      a = a * b;
      return a;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > operator / (Vector<VectorTag<AVX<double,SpT>,8> > const & a, Vector<VectorTag<AVX<double,SpT>,8> > const & b) {
      return _mm512_div_pd(a, b);
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > operator / (Vector<VectorTag<AVX<double,SpT>,8> > const & a, const double b) {
      return a / Vector<VectorTag<AVX<double,SpT>,8> >(b);
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > operator / (const double a, Vector<VectorTag<AVX<double,SpT>,8> > const & b) {
      return Vector<VectorTag<AVX<double,SpT>,8> >(a) / b;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > & operator /= (Vector<VectorTag<AVX<double,SpT>,8> > & a, Vector<VectorTag<AVX<double,SpT>,8> > const & b) {
      a = a / b;
      return a;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > & operator /= (Vector<VectorTag<AVX<double,SpT>,8> > & a, const double b) {
      a = a / b;
      return a;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<double,SpT>,8> > operator - (Vector<VectorTag<AVX<double,SpT>,8> > const & a) {
      return -1*a;
    }

    ///
    /// AVX512D Kokkos::complex<double>
    ///


    template<typename SpT>
    class Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > {
    public:
      using type = Vector<VectorTag<AVX<double,SpT>,4> >;
      using value_type = Kokkos::complex<double>;
      using real_type = double;

      enum : int { vector_length = 4 };

      union data_type {
        __m512d v;
        Kokkos::complex<double> d[4];

        data_type() { v = _mm512_setzero_pd(); }
        data_type(const data_type &b) { v = b.v; }
      };

    private:
      mutable data_type _data;

    public:
      inline Vector() { _data.v = _mm512_setzero_pd(); }
      inline Vector(const value_type val) {
        _data.v = _mm512_set_pd(val.imag(),val.real(),val.imag(),val.real(),
                                val.imag(),val.real(),val.imag(),val.real());
      }
      inline Vector(const real_type val) {
        _data.v = _mm512_set_pd(0,val,0,val,
                                0,val,0,val);
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
        _data.v = _mm512_load_pd((real_type*)p);
        return *this;
      }

      inline
      type& loadUnaligned(value_type const *p) {
        _data.v = _mm512_loadu_pd((real_type*)p);
        return *this;
      }

      inline
      void storeAligned(value_type *p) const {
        _mm512_store_pd((real_type*)p, _data.v);
      }

      inline
      void storeUnaligned(value_type *p) const {
        _mm512_storeu_pd((real_type*)p, _data.v);
      }

      inline
      value_type& operator[](int i) const {
        return _data.d[i];
      }

    };

    /// complex vector , complex vector

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator + (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & b) {
      return _mm512_add_pd(a, b);
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & operator += (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & b) {
      a = a + b;
      return a;
    }

    /// complex vector , complex

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator + (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & a, const Kokkos::complex<double> b) {
      return a + Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> >(b);
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator + (const Kokkos::complex<double> a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & b) {
      return Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> >(a) + b;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & operator += (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & a, const Kokkos::complex<double> b) {
      a = a + b;
      return a;
    }

    /// complex vector , real

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator + (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & a, const double b) {
      return a + Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> >(b);
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator + (const double a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & b) {
      return Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> >(a) + b;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & operator += (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & a, const double b) {
      a = a + b;
      return a;
    }

    /// unitary operator

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator ++ (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & a, int) {
      Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > a0 = a;
      a = a + 1.0;
      return a0;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & operator ++ (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & a) {
      a = a + 1.0;
      return a;
    }

    /// complex vector , complex vector

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator - (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & b) {
      return _mm512_sub_pd(a, b);
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & operator -= (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & b) {
      a = a - b;
      return a;
    }

    /// unitary operator

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator -- (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & a, int) {
      Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > a0 = a;
      a = a - 1.0;
      return a0;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & operator -- (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & a) {
      a = a - 1.0;
      return a;
    }

    /// complex vector , complex

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator - (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & a, const Kokkos::complex<double> b) {
      return a - Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> >(b);
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator - (const Kokkos::complex<double> a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & b) {
      return Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> >(a) - b;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & operator -= (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & a, const Kokkos::complex<double> b) {
      a = a - b;
      return a;
    }

    /// complex vector , real

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator - (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & a, const double b) {
      return a - Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> >(b);
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator - (const double a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & b) {
      return Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> >(a) - b;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & operator -= (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & a, const double b) {
      a = a - b;
      return a;
    }

    /// complex vector , complex vector

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator * (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & b) {
      const __m512d
        as = _mm512_permute_pd(a, 0x55),
        br = _mm512_permute_pd(b, 0x00),
        bi = _mm512_permute_pd(b, 0xff);

#if defined(__FMA__)
      return _mm512_fmaddsub_pd(a, br, _mm512_mul_pd(as, bi));
#else
      Kokkos::abort("Not yet implemented");
      return __m512d();
#endif
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & operator *= (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & b) {
      a = a * b;
      return a;
    }

    /// complex vector , complex

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator * (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & a, const Kokkos::complex<double> b) {
      return a * Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> >(b);
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator * (const Kokkos::complex<double> a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & b) {
      return Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> >(a) * b;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & operator *= (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & a, const Kokkos::complex<double> b) {
      a = a * b;
      return a;
    }

    /// complex vector , real

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator * (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & a, const double b) {
      return _mm512_mul_pd(a, _mm512_set1_pd(b));
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator * (const double a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & b) {
      return _mm512_mul_pd(_mm512_set1_pd(a), b);
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & operator *= (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & a, const double b) {
      a = a * b;
      return a;
    }

    /// complex vector , complex vector

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator / (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & b) {
      Kokkos::abort("Not yet implemented");
      return _mm512_div_pd(a, b);
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & operator /= (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & b) {
      Kokkos::abort("Not yet implemented");
      a = a / b;
      return a;
    }

    /// complex vector , complex

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator / (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & a, const Kokkos::complex<double> b) {
      Kokkos::abort("Not yet implemented");
      return a / Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> >(b);
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator / (const Kokkos::complex<double> a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & b) {
      Kokkos::abort("Not yet implemented");
      return Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> >(a) / b;
    }


    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & operator /= (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & a, const Kokkos::complex<double> b) {
      Kokkos::abort("Not yet implemented");
      a = a / b;
      return a;
    }

    /// complex vector , real

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator / (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & a, const double b) {
      Kokkos::abort("Not yet implemented");
      return a / Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> >(b);
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator / (const double a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & b) {
      Kokkos::abort("Not yet implemented");
      return Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> >(a) / b;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & operator /= (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > & a, const double b) {
      Kokkos::abort("Not yet implemented");
      a = a / b;
      return a;
    }

    template<typename SpT>
    inline
    static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator - (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & a) {
      return -1*a;
    }

  }
}
#endif
#endif
