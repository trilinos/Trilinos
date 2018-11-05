#ifndef __KOKKOSBATCHED_VECTOR_SIMD_HPP__
#define __KOKKOSBATCHED_VECTOR_SIMD_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <Kokkos_Complex.hpp>
#include <KokkosBatched_Vector.hpp>

#if defined(__CUDA_ARCH__) 
#undef  __KOKKOSBATCHED_ENABLE_AVX__
#else
// compiler bug with AVX in some architectures
#define __KOKKOSBATCHED_ENABLE_AVX__
#endif

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
      KOKKOS_INLINE_FUNCTION Vector(const ArgValueType &val) {
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
      KOKKOS_INLINE_FUNCTION Vector(const Vector<SIMD<ArgValueType>,vector_length> &b) {
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
      value_type& operator[](const int &i) const {
        return _data[i];
      }
    };

}}

#if defined(__KOKKOSBATCHED_ENABLE_AVX__)
#if defined(__AVX__) || defined(__AVX2__)
#include <immintrin.h>

namespace KokkosBatched {
  namespace Experimental {

    template<>
    class Vector<SIMD<double>,4> {
    public:
      using type = Vector<SIMD<double>,4>;
      using value_type = double;
      using mag_type = double;

      enum : int { vector_length = 4 };
      typedef __m256d data_type __attribute__ ((aligned(32)));

      KOKKOS_INLINE_FUNCTION
      static const char* label() { return "AVX256"; }

      template<typename,int>
      friend class Vector;

    private:
      mutable data_type _data;

    public:
      KOKKOS_INLINE_FUNCTION Vector() { _data = _mm256_setzero_pd(); }
      KOKKOS_INLINE_FUNCTION Vector(const value_type &val) { _data = _mm256_set1_pd(val); }
      KOKKOS_INLINE_FUNCTION Vector(const type &b) { _data = b._data; }
      KOKKOS_INLINE_FUNCTION Vector(const __m256d &val) { _data = val; }

      template<typename ArgValueType>
      KOKKOS_INLINE_FUNCTION Vector(const ArgValueType &val) {
        auto d = reinterpret_cast<value_type*>(&_data);
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
        for (int i=0;i<vector_length;++i)
          d[i] = val;
      }

      template<typename ArgValueType>
      KOKKOS_INLINE_FUNCTION Vector(const Vector<SIMD<ArgValueType>,vector_length> &b) {
	auto dd = reinterpret_cast<value_type*>(&_data);
	auto bb = reinterpret_cast<ArgValueType*>(&b._data);
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
	for (int i=0;i<vector_length;++i)
	  dd[i] = bb[i];
      }

      KOKKOS_INLINE_FUNCTION
      type& operator=(const __m256d &val) {
        _data = val;
        return *this;
      }

      KOKKOS_INLINE_FUNCTION
      operator __m256d() const {
        return _data;
      }

      KOKKOS_INLINE_FUNCTION
      type& loadAligned(const value_type *p) {
        _data = _mm256_load_pd(p);
        return *this;
      }

      KOKKOS_INLINE_FUNCTION
      type& loadUnaligned(const value_type *p) {
        _data = _mm256_loadu_pd(p);
        return *this;
      }

      KOKKOS_INLINE_FUNCTION
      void storeAligned(value_type *p) const {
        _mm256_store_pd(p, _data);
      }

      KOKKOS_INLINE_FUNCTION
      void storeUnaligned(value_type *p) const {
        _mm256_storeu_pd(p, _data);
      }

      KOKKOS_INLINE_FUNCTION
      value_type& operator[](const int &i) const {
        return reinterpret_cast<value_type*>(&_data)[i];
      }
    };

    template<>
    class Vector<SIMD<Kokkos::complex<double> >,2> {
    public:
      using type = Vector<SIMD<Kokkos::complex<double> >,2>;
      using value_type = Kokkos::complex<double>;
      using mag_type = double;

      static const int vector_length = 2;
      typedef __m256d data_type __attribute__ ((aligned(32)));

      KOKKOS_INLINE_FUNCTION
      static const char* label() { return "AVX256"; }
      
      template<typename,int>
      friend class Vector;

    private:
      mutable data_type _data;

    public:
      KOKKOS_INLINE_FUNCTION Vector() { _data = _mm256_setzero_pd(); }
      KOKKOS_INLINE_FUNCTION Vector(const value_type &val) { _data = _mm256_broadcast_pd((const __m128d *)&val);}
      KOKKOS_INLINE_FUNCTION Vector(const mag_type &val) { const value_type a(val); _data = _mm256_broadcast_pd((__m128d const *)&a); }
      KOKKOS_INLINE_FUNCTION Vector(const type &b) { _data = b._data; }
      KOKKOS_INLINE_FUNCTION Vector(const __m256d &val) { _data = val; }
      
//       template<typename ArgValueType>
//       KOKKOS_INLINE_FUNCTION Vector(const ArgValueType val) {
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
      KOKKOS_INLINE_FUNCTION Vector(const Vector<SIMD<ArgValueType>,vector_length> &b) {
	auto dd = reinterpret_cast<value_type*>(&_data);
	auto bb = reinterpret_cast<ArgValueType*>(&b._data);
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
	for (int i=0;i<vector_length;++i)
	  dd[i] = bb[i];
      }

      KOKKOS_INLINE_FUNCTION
      type& operator=(const __m256d &val) {
        _data = val;
        return *this;
      }

      KOKKOS_INLINE_FUNCTION
      operator __m256d() const {
        return _data;
      }

      KOKKOS_INLINE_FUNCTION
      type& loadAligned(const value_type *p) {
        _data = _mm256_load_pd((mag_type*)p);
        return *this;
      }

      KOKKOS_INLINE_FUNCTION
      type& loadUnaligned(const value_type *p) {
        _data = _mm256_loadu_pd((mag_type*)p);
        return *this;
      }

      KOKKOS_INLINE_FUNCTION
      void storeAligned(value_type *p) const {
        _mm256_store_pd((mag_type*)p, _data);
      }

      KOKKOS_INLINE_FUNCTION
      void storeUnaligned(value_type *p) const {
        _mm256_storeu_pd((mag_type*)p, _data);
      }

      KOKKOS_INLINE_FUNCTION
      value_type& operator[](const int &i) const {
        return reinterpret_cast<value_type*>(&_data)[i];
      }
    };
}}
#endif /* #if defined(__AVX__) || defined(__AVX2__) */

#if defined(__AVX512F__)
#include <immintrin.h>

namespace KokkosBatched {
  namespace Experimental {

    template<>
    class Vector<SIMD<double>,8> {
    public:
      using type = Vector<SIMD<double>,8>;
      using value_type = double;
      using mag_type = double;
      
      enum : int { vector_length = 8 };
      typedef __m512d data_type __attribute__ ((aligned(64)));

      template<typename,int>
      friend class Vector;

    private:
      mutable data_type _data;

    public:
      KOKKOS_INLINE_FUNCTION Vector() { _data = _mm512_setzero_pd(); }
      KOKKOS_INLINE_FUNCTION Vector(const value_type &val) { _data = _mm512_set1_pd(val); }
      KOKKOS_INLINE_FUNCTION Vector(const type &b) { _data = b._data; }
      KOKKOS_INLINE_FUNCTION Vector(const __m512d &val) { _data = val; }

      template<typename ArgValueType>
      KOKKOS_INLINE_FUNCTION Vector(const ArgValueType &val) {
        auto d = reinterpret_cast<value_type*>(&_data);
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
        for (int i=0;i<vector_length;++i)
          d[i] = val;
      }
      template<typename ArgValueType>
      KOKKOS_INLINE_FUNCTION Vector(const Vector<SIMD<ArgValueType>,vector_length> &b) {
	auto dd = reinterpret_cast<value_type*>(&_data);
	auto bb = reinterpret_cast<ArgValueType*>(&b._data);
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
	for (int i=0;i<vector_length;++i)
	  dd[i] = bb[i];
      }

      KOKKOS_INLINE_FUNCTION
      type& operator=(const __m512d &val) {
        _data = val;
        return *this;
      }

      KOKKOS_INLINE_FUNCTION
      operator __m512d() const {
        return _data;
      }

      KOKKOS_INLINE_FUNCTION
      type& loadAligned(const value_type *p) {
        _data = _mm512_load_pd(p);
        return *this;
      }

      KOKKOS_INLINE_FUNCTION
      type& loadUnaligned(const value_type *p) {
        _data = _mm512_loadu_pd(p);
        return *this;
      }

      KOKKOS_INLINE_FUNCTION
      void storeAligned(value_type *p) const {
        _mm512_store_pd(p, _data);
      }

      KOKKOS_INLINE_FUNCTION
      void storeUnaligned(value_type *p) const {
        _mm512_storeu_pd(p, _data);
      }

      KOKKOS_INLINE_FUNCTION
      value_type& operator[](const int &i) const {
        return reinterpret_cast<value_type*>(&_data)[i];
      }
    };

    template<>
    class Vector<SIMD<Kokkos::complex<double> >,4> {
    public:
      using type = Vector<SIMD<Kokkos::complex<double> >,4>;
      using value_type = Kokkos::complex<double>;
      using mag_type = double;

      enum : int { vector_length = 4 };
      typedef __m512d data_type __attribute__ ((aligned(64)));

      KOKKOS_INLINE_FUNCTION
      static const char* label() { return "AVX512"; }

      template<typename,int>
      friend class Vector;

    private:
      mutable data_type _data;

    public:
      KOKKOS_INLINE_FUNCTION Vector() { _data = _mm512_setzero_pd(); }
      KOKKOS_INLINE_FUNCTION Vector(const value_type &val) {
        _data = _mm512_mask_broadcast_f64x4(_mm512_set1_pd(val.imag()), 0x55, _mm256_set1_pd(val.real()));
      }
      KOKKOS_INLINE_FUNCTION Vector(const mag_type &val) {
        _data = _mm512_mask_broadcast_f64x4(_mm512_setzero_pd(), 0x55, _mm256_set1_pd(val));
      }
      KOKKOS_INLINE_FUNCTION Vector(const type &b) { _data = b._data; }
      KOKKOS_INLINE_FUNCTION Vector(const __m512d &val) { _data = val; }

      template<typename ArgValueType>
      KOKKOS_INLINE_FUNCTION Vector(const ArgValueType &val) {
        auto d = reinterpret_cast<value_type*>(&_data);
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
        for (int i=0;i<vector_length;++i)
          d[i] = val;
      }
      template<typename ArgValueType>
      KOKKOS_INLINE_FUNCTION Vector(const Vector<SIMD<ArgValueType>,vector_length> &b) {
	auto dd = reinterpret_cast<value_type*>(&_data);
	auto bb = reinterpret_cast<value_type*>(&b._data);
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
	for (int i=0;i<vector_length;++i)
	  dd[i] = bb[i];
      }

      KOKKOS_INLINE_FUNCTION
      type& operator=(const __m512d &val) {
        _data = val;
        return *this;
      }

      KOKKOS_INLINE_FUNCTION
      operator __m512d() const {
        return _data;
      }

      KOKKOS_INLINE_FUNCTION
      type& loadAligned(const value_type *p) {
        _data = _mm512_load_pd((mag_type*)p);
        return *this;
      }

      KOKKOS_INLINE_FUNCTION
      type& loadUnaligned(const value_type *p) {
        _data = _mm512_loadu_pd((mag_type*)p);
        return *this;
      }

      KOKKOS_INLINE_FUNCTION
      void storeAligned(value_type *p) const {
        _mm512_store_pd((mag_type*)p, _data);
      }

      KOKKOS_INLINE_FUNCTION
      void storeUnaligned(value_type *p) const {
        _mm512_storeu_pd((mag_type*)p, _data);
      }

      KOKKOS_INLINE_FUNCTION
      value_type& operator[](const int &i) const {
        return reinterpret_cast<value_type*>(&_data)[i];
      }
    };
}}
#endif /* #if defined(__AVX512F__) */
#endif /* #if defined(__KOKKOSBATCHED_ENABLE_AVX__) */

#include "KokkosBatched_Vector_SIMD_Arith.hpp"
#include "KokkosBatched_Vector_SIMD_Logical.hpp"
#include "KokkosBatched_Vector_SIMD_Relation.hpp"
#include "KokkosBatched_Vector_SIMD_Math.hpp"
#include "KokkosBatched_Vector_SIMD_Misc.hpp"
#include "KokkosBatched_Vector_SIMD_View.hpp"

#endif
