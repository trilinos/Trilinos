#ifndef __KOKKOSKERNELS_VECTOR_AVX256D_HPP__
#define __KOKKOSKERNELS_VECTOR_AVX256D_HPP__

/// This follows the work "vectorlib" -- GNU public license -- written by Agner Fog.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#if defined(__AVX__) || defined(__AVX2__) 

#include <immintrin.h>

namespace KokkosKernels {

  ///
  /// AVX256D double
  ///

  template<typename SpT>
  class Vector<VectorTag<AVX<double,SpT>,4> > {
  public:
    using type = Vector<VectorTag<AVX<double,SpT>,4> >;
    using value_type = double;
    using real_type = double;

    enum : int { vector_length = 4 };
 
    union data_type {
      __m256d v;
      double d[4];
    };

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

    template <int i0, int i1, int i2, int i3, int i4, int i5, int i6, int i7>
    static inline __m256 constant8f() {
      static const union {
        int  i[8];
      __m256 ymm;
      } u = {{i0,i1,i2,i3,i4,i5,i6,i7}};
      return u.ymm;
    }
    
  };

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > operator + (Vector<VectorTag<AVX<double,SpT>,4> > const & a, Vector<VectorTag<AVX<double,SpT>,4> > const & b) {
    return _mm256_add_pd(a, b);
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > operator + (Vector<VectorTag<AVX<double,SpT>,4> > const & a, const double b) {
    return a + Vector<VectorTag<AVX<double,SpT>,4> >(b);
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > operator + (const double a, Vector<VectorTag<AVX<double,SpT>,4> > const & b) {
    return Vector<VectorTag<AVX<double,SpT>,4> >(a) + b;
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > & operator += (Vector<VectorTag<AVX<double,SpT>,4> > & a, Vector<VectorTag<AVX<double,SpT>,4> > const & b) {
    a = a + b;
    return a;
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > & operator += (Vector<VectorTag<AVX<double,SpT>,4> > & a, const double b) {
    a = a + b;
    return a;
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > operator ++ (Vector<VectorTag<AVX<double,SpT>,4> > & a, int) {
    Vector<VectorTag<AVX<double,SpT>,4> > a0 = a;
    a = a + 1.0;
    return a0;
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > & operator ++ (Vector<VectorTag<AVX<double,SpT>,4> > & a) {
    a = a + 1.0;
    return a;
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > operator - (Vector<VectorTag<AVX<double,SpT>,4> > const & a, Vector<VectorTag<AVX<double,SpT>,4> > const & b) {
    return _mm256_sub_pd(a, b);
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > operator - (Vector<VectorTag<AVX<double,SpT>,4> > const & a, const double b) {
    return a - Vector<VectorTag<AVX<double,SpT>,4> >(b);
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > operator - (const double a, Vector<VectorTag<AVX<double,SpT>,4> > const & b) {
    return Vector<VectorTag<AVX<double,SpT>,4> >(a) - b;
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > operator - (Vector<VectorTag<AVX<double,SpT>,4> > const & a) {
    return _mm256_xor_pd(a, _mm256_castps_pd(Vector<VectorTag<AVX<double,SpT>,4> >::constant8f<
                                             0,(int)0x80000000,
                                             0,(int)0x80000000,
                                             0,(int)0x80000000,
                                             0,(int)0x80000000> ()));
  }
  
  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > & operator -= (Vector<VectorTag<AVX<double,SpT>,4> > & a, Vector<VectorTag<AVX<double,SpT>,4> > const & b) {
    a = a - b;
    return a;
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > & operator -= (Vector<VectorTag<AVX<double,SpT>,4> > & a, const double b) {
    a = a - b;
    return a;
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > operator -- (Vector<VectorTag<AVX<double,SpT>,4> > & a, int) {
    Vector<VectorTag<AVX<double,SpT>,4> > a0 = a;
    a = a - 1.0;
    return a0;
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > & operator -- (Vector<VectorTag<AVX<double,SpT>,4> > & a) {
    a = a - 1.0;
    return a;
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > operator * (Vector<VectorTag<AVX<double,SpT>,4> > const & a, Vector<VectorTag<AVX<double,SpT>,4> > const & b) {
    return _mm256_mul_pd(a, b);
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > operator * (Vector<VectorTag<AVX<double,SpT>,4> > const & a, const double b) {
    return a * Vector<VectorTag<AVX<double,SpT>,4> >(b);
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > operator * (const double a, Vector<VectorTag<AVX<double,SpT>,4> > const & b) {
    return Vector<VectorTag<AVX<double,SpT>,4> >(a) * b;
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > & operator *= (Vector<VectorTag<AVX<double,SpT>,4> > & a, Vector<VectorTag<AVX<double,SpT>,4> > const & b) {
    a = a * b;
    return a;
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > & operator *= (Vector<VectorTag<AVX<double,SpT>,4> > & a, const double b) {
    a = a * b;
    return a;
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > operator / (Vector<VectorTag<AVX<double,SpT>,4> > const & a, Vector<VectorTag<AVX<double,SpT>,4> > const & b) {
    return _mm256_div_pd(a, b);
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > operator / (Vector<VectorTag<AVX<double,SpT>,4> > const & a, const double b) {
    return a / Vector<VectorTag<AVX<double,SpT>,4> >(b);
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > operator / (const double a, Vector<VectorTag<AVX<double,SpT>,4> > const & b) {
    return Vector<VectorTag<AVX<double,SpT>,4> >(a) / b;
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > & operator /= (Vector<VectorTag<AVX<double,SpT>,4> > & a, Vector<VectorTag<AVX<double,SpT>,4> > const & b) {
    a = a / b;
    return a;
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<double,SpT>,4> > & operator /= (Vector<VectorTag<AVX<double,SpT>,4> > & a, const double b) {
    a = a / b;
    return a;
  }

  ///
  /// AVX256D Kokkos::complex<double>
  //

  template<typename SpT>
  class Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > {
  public:
    using type = Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> >;
    using value_type = Kokkos::complex<double>;
    using real_type = double;

    enum : int { vector_length = 2 };
 
    union data_type {
      __m256d v;
      Kokkos::complex<double> d[2];

      data_type() { v = _mm256_setzero_pd(); }
      data_type(const data_type &b) { v = b.v; }
    };

  private:
    mutable data_type _data;
    
  public:
    inline Vector() { _data.v = _mm256_setzero_pd(); }
    inline Vector(const value_type val) { _data.v = _mm256_set_pd(val.imag(),val.real(),val.imag(),val.real()); }
    inline Vector(const real_type val) { _data.v = _mm256_set_pd(0, val, 0, val); }
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
      _data.v = _mm256_load_pd((real_type*)p);
      return *this;
    }

    inline 
    type& loadUnaligned(value_type const *p) {
      _data.v = _mm256_loadu_pd((real_type*)p);
      return *this;
    }
    
    inline 
    void storeAligned(value_type *p) const {
      _mm256_store_pd((real_type*)p, _data.v);
    }
    
    inline 
    void storeUnaligned(value_type *p) const {
      _mm256_storeu_pd((real_type*)p, _data.v);
    }

    inline 
    value_type& operator[](int i) const {
      return _data.d[i];
    }

    template <int i0, int i1, int i2, int i3, int i4, int i5, int i6, int i7>
    static inline __m256 constant8f() {
      static const union {
        int  i[8];
      __m256 ymm;
      } u = {{i0,i1,i2,i3,i4,i5,i6,i7}};
      return u.ymm;
    }
    
  };

  /// complex vector , complex vector

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator + (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & b) {
    return _mm256_add_pd(a, b);
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & operator += (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & b) {
    a = a + b;
    return a;
  }

  /// complex vector , complex

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator + (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & a, const Kokkos::complex<double> b) {
    return a + Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> >(b);
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator + (const Kokkos::complex<double> a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & b) {
    return Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> >(a) + b;
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & operator += (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & a, const Kokkos::complex<double> b) {
    a = a + b;
    return a;
  }

  /// complex vector , real

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator + (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & a, const double b) {
    return a + Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> >(b);
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator + (const double a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & b) {
    return Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> >(a) + b;
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & operator += (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & a, const double b) {
    a = a + b;
    return a;
  }

  /// unitary operator

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator ++ (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & a, int) {
    Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > a0 = a;
    a = a + 1.0;
    return a0;
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & operator ++ (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & a) {
    a = a + 1.0;
    return a;
  }

  /// complex vector , complex vector

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator - (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & b) {
    return _mm256_sub_pd(a, b);
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & operator -= (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & b) {
    a = a - b;
    return a;
  }

  /// complex vector , complex

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator - (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & a, const Kokkos::complex<double> b) {
    return a - Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> >(b);
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator - (const Kokkos::complex<double> a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & b) {
    return Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> >(a) - b;
  }
  
  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & operator -= (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & a, const Kokkos::complex<double> b) {
    a = a - b;
    return a;
  }

  /// complex vector , real

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator - (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & a, const double b) {
    return a - Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> >(b);
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator - (const double a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & b) {
    return Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> >(a) - b;
  }
  
  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & operator -= (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & a, const double b) {
    a = a - b;
    return a;
  }

  /// unitary operator

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator - (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & a) {
    return _mm256_xor_pd(a, _mm256_castps_pd(Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> >::constant8f<
                                             0,(int)0x80000000,
                                             0,(int)0x80000000,
                                             0,(int)0x80000000,
                                             0,(int)0x80000000> ()));
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator -- (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & a, int) {
    Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > a0 = a;
    a = a - 1.0;
    return a0;
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & operator -- (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & a) {
    a = a - 1.0;
    return a;
  }

  /// complex vector , complex vector

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator * (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & b) {
    const __m256d 
      as = _mm256_permute_pd(a, 0x5),
      br = _mm256_permute_pd(b, 0x0),
      bi = _mm256_permute_pd(b, 0xf);

    return _mm256_fmaddsub_pd(a, br, _mm256_mul_pd(as, bi));
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & operator *= (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & b) {
    a = a * b;
    return a;
  }

  /// complex vector , complex 

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator * (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & a, const Kokkos::complex<double> b) {
    return a * Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> >(b);
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator * (const Kokkos::complex<double> a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & b) {
    return Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> >(a) * b;
  }
  
  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & operator *= (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & a, const Kokkos::complex<double> b) {
    a = a * b;
    return a;
  }

  /// complex vector , real

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator * (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & a, const double b) {
    return _mm256_mul_pd(a, _mm256_set1_pd(b)); 
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator * (const double a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & b) {
    return _mm256_mul_pd(_mm256_set1_pd(a), b); 
  }
  
  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & operator *= (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & a, const double b) {
    a = a * b;
    return a;
  }

  /// complex vector , complex vector

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator / (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & b) {
    Kokkos::abort("Not yet implemented");
    return _mm256_div_pd(a, b);
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & operator /= (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & b) {
    Kokkos::abort("Not yet implemented");
    a = a / b;
    return a;
  }

  /// complex vector , complex
  
  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator / (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & a, const Kokkos::complex<double> b) {
    Kokkos::abort("Not yet implemented");
    return a / Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> >(b);
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator / (const Kokkos::complex<double> a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & b) {
    Kokkos::abort("Not yet implemented");
    return Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> >(a) / b;
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & operator /= (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & a, const Kokkos::complex<double> b) {
    Kokkos::abort("Not yet implemented");
    a = a / b;
    return a;
  }

  /// complex vector , real
  
  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator / (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & a, const double b) {
    Kokkos::abort("Not yet implemented");
    return a / Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> >(b);
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > operator / (const double a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > const & b) {
    Kokkos::abort("Not yet implemented");
    return Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> >(a) / b;
  }

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & operator /= (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,2> > & a, const double b) {
    Kokkos::abort("Not yet implemented");
    a = a / b;
    return a;
  }

}

#endif
#endif
