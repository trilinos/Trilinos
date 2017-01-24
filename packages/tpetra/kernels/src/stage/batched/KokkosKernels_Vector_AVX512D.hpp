/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels: Linear Algebra and Graph Kernels
//                 Copyright 2016 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/





#ifndef __KOKKOSKERNELS_VECTOR_AVX512D_HPP__
#define __KOKKOSKERNELS_VECTOR_AVX512D_HPP__

/// This follows the work "vectorlib" -- GNU public license -- written by Agner Fog.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#if defined(__AVX__) || defined(__AVX2__)

#include <immintrin.h>

namespace KokkosKernels {

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

    template <int i00, int i01, int i02, int i03, int i04, int i05, int i06, int i07,
              int i10, int i11, int i12, int i13, int i14, int i15, int i16, int i17>
    static inline __m512 constant16f() {
      static const union {
        int  i[16];
      __m512 zmm;
      } u = {{i00,i01,i02,i03,i04,i05,i06,i07,
              i10,i11,i12,i13,i14,i15,i16,i17}};
      return u.zmm;
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
  static Vector<VectorTag<AVX<double,SpT>,8> > operator - (Vector<VectorTag<AVX<double,SpT>,8> > const & a) {
    return _mm512_xor_pd(a, _mm512_castps_pd(Vector<VectorTag<AVX<double,SpT>,8> >::constant16f<
                                             0,(int)0x80000000,
                                             0,(int)0x80000000,
                                             0,(int)0x80000000,
                                             0,(int)0x80000000,
                                             0,(int)0x80000000,
                                             0,(int)0x80000000,
                                             0,(int)0x80000000,
                                             0,(int)0x80000000> ()));
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

    template <int i00, int i01, int i02, int i03, int i04, int i05, int i06, int i07,
              int i10, int i11, int i12, int i13, int i14, int i15, int i16, int i17>
    static inline __m512 constant16f() {
      static const union {
        int  i[16];
      __m512 zmm;
      } u = {{i00,i01,i02,i03,i04,i05,i06,i07,
              i10,i11,i12,i13,i14,i15,i16,i17}};
      return u.zmm;
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

  /// unitary operator

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator - (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & a) {
    return _mm512_xor_pd(a, _mm512_castps_pd(Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> >::constant16f<
                                             0,(int)0x80000000,
                                             0,(int)0x80000000,
                                             0,(int)0x80000000,
                                             0,(int)0x80000000,
                                             0,(int)0x80000000,
                                             0,(int)0x80000000,
                                             0,(int)0x80000000,
                                             0,(int)0x80000000> ()));
  }

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

  /// complex vector , complex vector

  template<typename SpT>
  inline 
  static Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > operator * (Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & a, Vector<VectorTag<AVX<Kokkos::complex<double>,SpT>,4> > const & b) {
    const __m512d
      as = _mm512_permute_pd(a, 0x55),
      br = _mm512_permute_pd(b, 0x00),
      bi = _mm512_permute_pd(b, 0xff);

    return _mm512_fmaddsub_pd(a, br, _mm512_mul_pd(as, bi));
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


}

#endif
#endif
