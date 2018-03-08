#ifndef __KOKKOSBATCHED_VECTOR_SIMD_HPP__
#define __KOKKOSBATCHED_VECTOR_SIMD_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <Kokkos_Complex.hpp>
#include <KokkosBatched_Vector.hpp>

namespace KokkosBatched {
  namespace Experimental {

    template<typename T, int v = 0>
    struct TypeTraits {
      typedef T thread_private_type;
      typedef T team_private_type;
    };

    template<typename T, int l, int v>
    struct TypeTraits<Vector<SIMD<T>,l>, v> {
      typedef typename std::conditional<std::is_same<Kokkos::Impl::ActiveExecutionMemorySpace,Kokkos::HostSpace>::value, 
                                        Vector<SIMD<T>,l>, T>::type thread_private_type;
      typedef typename std::conditional<std::is_same<Kokkos::Impl::ActiveExecutionMemorySpace,Kokkos::HostSpace>::value, 
                                        Vector<SIMD<T>,l>, Vector<SIMD<T>,(l/v)+(l%v>0)> > team_private_type;      
    };
  
    template<typename T, int l>
    class Vector<SIMD<T>,l> {
    public:
      using type = Vector<SIMD<T>,l>;
      using value_type = T;
      using mag_type = typename Kokkos::Details::ArithTraits<T>::mag_type;

      enum : int { vector_length = l };
      typedef value_type data_type[vector_length];

      KOKKOS_INLINE_FUNCTION
      static const char* label() { return "SIMD"; }

      template<typename,int>
      friend class Vector;

    private:
      mutable data_type _data;

    public:
      KOKKOS_INLINE_FUNCTION Vector() {
        //static_assert(std::is_same<Kokkos::Impl::ActiveExecutionMemorySpace,Kokkos::HostSpace>::value,
        //              "Vector SIMD should not be instanciated in CudaSpace");
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
        for (int i=0;i<vector_length;++i)
          _data[i] = 0;
      }
      template<typename ArgValueType>
      KOKKOS_INLINE_FUNCTION Vector(const ArgValueType val) {
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
        for (int i=0;i<vector_length;++i)
          _data[i] = val;
      }
      template<typename ArgValueType>
      KOKKOS_INLINE_FUNCTION Vector(const Vector<SIMD<ArgValueType>,vector_length> b) {
        static_assert(std::is_convertible<value_type,ArgValueType>::value, "input type is not convertible");
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
        for (int i=0;i<vector_length;++i)
          _data[i] = b[i];
      }

      KOKKOS_INLINE_FUNCTION
      type& loadAligned(const value_type *p) {
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
        for (int i=0;i<vector_length;++i)
          _data[i] = p[i];
        return *this;
      }

      KOKKOS_INLINE_FUNCTION
      type& loadUnaligned(const value_type *p) {
        return loadAligned(p);
      }
      
      KOKKOS_INLINE_FUNCTION
      void storeAligned(value_type *p) const {
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
        for (int i=0;i<vector_length;++i)
          p[i] = _data[i];
      }

      KOKKOS_INLINE_FUNCTION
      void storeUnaligned(value_type *p) const {
        storeAligned(p);
      }

      KOKKOS_INLINE_FUNCTION
      value_type& operator[](const int i) const {
        return _data[i];
      }
    };

#ifndef __CUDA_ARCH__
#if defined(__AVX__) || defined(__AVX2__)
#include <immintrin.h>

    template<>
    class Vector<SIMD<double>,4> {
    public:
      using type = Vector<SIMD<double>,4>;
      using value_type = double;
      using mag_type = double;

      enum : int { vector_length = 4 };

      union data_type {
        __m256d v;
        double d[4];
      };

      inline
      static const char* label() { return "AVX256"; }

      template<typename,int>
      friend class Vector;

    private:
      mutable data_type _data;

    public:
      inline Vector() { _data.v = _mm256_setzero_pd(); }
      inline Vector(const value_type val) { _data.v = _mm256_set1_pd(val); }
      inline Vector(const type &b) { _data.v = b._data.v; }
      inline Vector(const __m256d &val) { _data.v = val; }

      template<typename ArgValueType>
      inline Vector(const ArgValueType val) {
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
        for (int i=0;i<vector_length;++i)
          _data.d[i] = val;
      }
      template<typename ArgValueType>
      inline Vector(const Vector<SIMD<ArgValueType>,vector_length> b) {
        static_assert(std::is_convertible<value_type,ArgValueType>::value, "input type is not convertible");
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
        for (int i=0;i<vector_length;++i)
          _data.d[i] = b[i];
      }

      inline
      type& operator=(const __m256d &val) {
        _data.v = val;
        return *this;
      }

      inline
      operator __m256d() const {
        return _data.v;
      }

      inline
      type& loadAligned(const value_type *p) {
        _data.v = _mm256_load_pd(p);
        return *this;
      }

      inline
      type& loadUnaligned(const value_type *p) {
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

      KOKKOS_INLINE_FUNCTION
      value_type& operator[](const int i) const {
        return _data.d[i];
      }
    };

    template<>
    class Vector<SIMD<Kokkos::complex<double> >,2> {
    public:
      using type = Vector<SIMD<Kokkos::complex<double> >,2>;
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

      template<typename,int>
      friend class Vector;

    private:
      mutable data_type _data;

    public:
      inline Vector() { _data.v = _mm256_setzero_pd(); }
      inline Vector(const value_type val) { _data.v = _mm256_broadcast_pd((const __m128d *)&val);}
      inline Vector(const mag_type val) { const value_type a(val); _data.v = _mm256_broadcast_pd((__m128d const *)&a); }
      inline Vector(const type &b) { _data.v = b._data.v; }
      inline Vector(const __m256d &val) { _data.v = val; }

//       template<typename ArgValueType>
//       inline Vector(const ArgValueType val) {
// #if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
// #pragma ivdep
// #endif
// #if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
// #pragma vector always
// #endif
//         for (int i=0;i<vector_length;++i)
//           _data.d[i] = value_type(val);
//       }
      template<typename ArgValueType>
      inline Vector(const Vector<SIMD<ArgValueType>,vector_length> b) {
        static_assert(std::is_convertible<value_type,ArgValueType>::value, "input type is not convertible");
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
        for (int i=0;i<vector_length;++i)
          _data.d[i] = b[i];
      }

      inline
      type& operator=(const __m256d &val) {
        _data.v = val;
        return *this;
      }

      inline
      operator __m256d() const {
        return _data.v;
      }

      inline
      type& loadAligned(const value_type *p) {
        _data.v = _mm256_load_pd((mag_type*)p);
        return *this;
      }

      inline
      type& loadUnaligned(const value_type *p) {
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

      KOKKOS_INLINE_FUNCTION
      value_type& operator[](int i) const {
        return _data.d[i];
      }
    };
#endif

#if defined(__AVX512F__)
#include <immintrin.h>
    template<>
    class Vector<SIMD<double>,8> {
    public:
      using type = Vector<SIMD<double>,8>;
      using value_type = double;
      using mag_type = double;
      
      enum : int { vector_length = 8 };

      union data_type {
        __m512d v;
        double d[8];
      };

      KOKKOS_INLINE_FUNCTION
      static const char* label() { return "AVX512"; }

      template<typename,int>
      friend class Vector;

    private:
      mutable data_type _data;

    public:
      inline Vector() { _data.v = _mm512_setzero_pd(); }
      inline Vector(const value_type val) { _data.v = _mm512_set1_pd(val); }
      inline Vector(const type &b) { _data.v = b._data.v; }
      inline Vector(const __m512d &val) { _data.v = val; }

      template<typename ArgValueType>
      inline Vector(const ArgValueType val) {
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
        for (int i=0;i<vector_length;++i)
          _data.d[i] = val;
      }
      template<typename ArgValueType>
      inline Vector(const Vector<SIMD<ArgValueType>,vector_length> b) {
        static_assert(std::is_convertible<value_type,ArgValueType>::value, "input type is not convertible");
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
        for (int i=0;i<vector_length;++i)
          _data.d[i] = b[i];
      }

      inline
      type& operator=(const __m512d &val) {
        _data.v = val;
        return *this;
      }

      inline
      operator __m512d() const {
        return _data.v;
      }

      inline
      type& loadAligned(const value_type *p) {
        _data.v = _mm512_load_pd(p);
        return *this;
      }

      inline
      type& loadUnaligned(const value_type *p) {
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

      KOKKOS_INLINE_FUNCTION
      value_type& operator[](const int i) const {
        return _data.d[i];
      }
    };

    template<>
    class Vector<SIMD<Kokkos::complex<double> >,4> {
    public:
      using type = Vector<SIMD<Kokkos::complex<double> >,4>;
      using value_type = Kokkos::complex<double>;
      using mag_type = double;

      enum : int { vector_length = 4 };

      union data_type {
        __m512d v;
        Kokkos::complex<double> d[4];

        data_type() { v = _mm512_setzero_pd(); }
        data_type(const data_type &b) { v = b.v; }
      };

      KOKKOS_INLINE_FUNCTION
      static const char* label() { return "AVX512"; }

      template<typename,int>
      friend class Vector;

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
      inline Vector(const __m512d &val) { _data.v = val; }

      template<typename ArgValueType>
      inline Vector(const ArgValueType val) {
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
        for (int i=0;i<vector_length;++i)
          _data.d[i] = val;
      }
      template<typename ArgValueType>
      inline Vector(const Vector<SIMD<ArgValueType>,vector_length> b) {
        static_assert(std::is_convertible<value_type,ArgValueType>::value, "input type is not convertible");
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
        for (int i=0;i<vector_length;++i)
          _data.d[i] = b[i];
      }

      inline
      type& operator=(const __m512d &val) {
        _data.v = val;
        return *this;
      }

      inline
      operator __m512d() const {
        return _data.v;
      }

      inline
      type& loadAligned(const value_type *p) {
        _data.v = _mm512_load_pd((mag_type*)p);
        return *this;
      }

      inline
      type& loadUnaligned(const value_type *p) {
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

      KOKKOS_INLINE_FUNCTION
      value_type& operator[](const int i) const {
        return _data.d[i];
      }
    };
#endif
#endif
  }
}

#include "KokkosBatched_Vector_SIMD_Arith.hpp"
#include "KokkosBatched_Vector_SIMD_Logical.hpp"
#include "KokkosBatched_Vector_SIMD_Relation.hpp"
#include "KokkosBatched_Vector_SIMD_Math.hpp"
#include "KokkosBatched_Vector_SIMD_Misc.hpp"
#include "KokkosBatched_Vector_SIMD_View.hpp"

#endif
