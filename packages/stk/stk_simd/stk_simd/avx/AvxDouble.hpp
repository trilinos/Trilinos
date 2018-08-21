// Copyright 2013 Sandia Corporation, Albuquerque, NM.

// IWYU pragma: private, include <stk_simd/Simd.hpp>

namespace stk {
namespace simd {

struct Double {

  STK_MATH_FORCE_INLINE Double() {}


  template <typename T>
  STK_MATH_FORCE_INLINE Double(const T x, typename std::enable_if<std::is_convertible<T,double>::value, void*>::type=0)
    : _data(_mm256_set1_pd(double(x))) {
  }

  STK_MATH_FORCE_INLINE Double(const __m256d& x) 
    : _data(x) {
  }

  STK_MATH_FORCE_INLINE Double(const Double& x) 
    : _data(x._data) {
  }

  STK_MATH_FORCE_INLINE Double(const volatile Double& x)
    : _data(x._data) {
  }

  STK_MATH_FORCE_INLINE Double& operator= (const Double& x) {
    _data = x._data;
    return *this;
  }

  template <typename T>
  STK_MATH_FORCE_INLINE typename std::enable_if<std::is_convertible<T,double>::value, Double&>::type operator= (const T x) {
    _data = _mm256_set1_pd(double(x));
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator+= (const Double& a) {
    _data = _mm256_add_pd(_data,a._data);
    return *this;
  }

  STK_MATH_FORCE_INLINE void operator+= (const volatile Double& a) volatile {
    _data = _mm256_add_pd(_data,a._data);
  }

  STK_MATH_FORCE_INLINE Double& operator-= (const Double& a) {
    _data = _mm256_sub_pd(_data,a._data);
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator*= (const Double& a) {
    _data = _mm256_mul_pd(_data,a._data);
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator/= (const Double& a) {
    _data = _mm256_div_pd(_data,a._data);
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator+= (const double a) {
    _data = _mm256_add_pd(_data,_mm256_set1_pd(a));
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator-= (const double a) {
    _data = _mm256_sub_pd(_data,_mm256_set1_pd(a));
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator*= (const double a) {
    _data = _mm256_mul_pd(_data,_mm256_set1_pd(a));
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator/= (const double a) {
    _data = _mm256_div_pd(_data,_mm256_set1_pd(a));
    return *this;
  }

  STK_MATH_FORCE_INLINE Double operator-() const {
    return Double( _mm256_sub_pd(Double(0.0)._data,_data) );
  }

  STK_MATH_FORCE_INLINE double& operator[](int i) {return (reinterpret_cast<double*>(&_data))[i];}
  STK_MATH_FORCE_INLINE const double& operator[](int i) const {return (reinterpret_cast<const double*>(&_data))[i];}
    
  STK_MATH_FORCE_INLINE int64_t& Int(int i) {return (reinterpret_cast<int64_t*>(&_data))[i];}
  STK_MATH_FORCE_INLINE const int64_t& Int(int i) const {return (reinterpret_cast<const int64_t*>(&_data))[i];}

  STK_MATH_FORCE_INLINE uint64_t& UInt(int i) {return (reinterpret_cast<uint64_t*>(&_data))[i];}
  STK_MATH_FORCE_INLINE const uint64_t& UInt(int i) const {return (reinterpret_cast<const uint64_t*>(&_data))[i];}

  __m256d _data; // the "_" means you should try not to use this directly
  // it is made public to avoid function call overhead 
  // and/or so the compiler doesn't have to use up one of
  // inlining depths (usually max inlining depth ~5)

};

} // namespace simd
} // namespace stk

