// Copyright 2013 Sandia Corporation, Albuquerque, NM.
// IWYU pragma: private, include <stk_simd/Simd.hpp>

#include <array>
#include "./FwSizes.hpp"

namespace stk {
namespace simd {

struct Float {
  
  STK_MATH_FORCE_INLINE Float() {}

  template <typename T>
  STK_MATH_FORCE_INLINE Float(const T x, typename std::enable_if<std::is_convertible<T,float>::value, void*>::type=0) {
    for (int i=0; i < nfloats; ++i) _data[i] = float(x);
  }

  STK_MATH_FORCE_INLINE Float(const Float& x) {
    for (int i=0; i < nfloats; ++i) _data[i] = x._data[i];
  }

  STK_MATH_FORCE_INLINE Float& operator= (const Float& x) {
    for (int i=0; i < nfloats; ++i) _data[i] = x._data[i];
    return *this;
  }

  template <typename T>
  STK_MATH_FORCE_INLINE typename std::enable_if<std::is_convertible<T,float>::value, Float&>::type operator= (const float x) {
    for (int i=0; i < nfloats; ++i) _data[i] = float(x);
    return *this;
  }

  STK_MATH_FORCE_INLINE Float& operator+= (const Float& a) {
    for (int i=0; i < nfloats; ++i) _data[i] += a[i];
    return *this;
  }
  
  STK_MATH_FORCE_INLINE Float& operator-= (const Float& a) {
    for (int i=0; i < nfloats; ++i) _data[i] -= a[i];
    return *this;
  }

  STK_MATH_FORCE_INLINE Float& operator*= (const Float& a) {
    for (int i=0; i < nfloats; ++i) _data[i] *= a[i];
    return *this;
  }

  STK_MATH_FORCE_INLINE Float& operator/= (const Float& a) {
    for (int i=0; i < nfloats; ++i) _data[i] /= a[i];
    return *this;
  }

  STK_MATH_FORCE_INLINE Float& operator+= (const float a) {
    for (int i=0; i < nfloats; ++i) _data[i] += a;
    return *this;
  }

  STK_MATH_FORCE_INLINE Float& operator-= (const float a) {
    for (int i=0; i < nfloats; ++i) _data[i] -= a;
    return *this;
  }

  STK_MATH_FORCE_INLINE Float& operator*= (const float a) {
    for (int i=0; i < nfloats; ++i) _data[i] *= a;
    return *this;
  }

  STK_MATH_FORCE_INLINE Float& operator/= (const float a) {
    for (int i=0; i < nfloats; ++i) _data[i] /= a;
    return *this;
  }

  STK_MATH_FORCE_INLINE Float operator-() const {
    Float tmp;
    for (int i=0; i < nfloats; ++i) tmp[i] = -_data[i];
    return tmp;
  }

  STK_MATH_FORCE_INLINE float& operator[](int i) {return _data[i];}
  STK_MATH_FORCE_INLINE const float& operator[](int i) const {return _data[i];}
    
  STK_MATH_FORCE_INLINE int32_t& Int(int i) {return (reinterpret_cast<int32_t*>(&_data))[i];}
  STK_MATH_FORCE_INLINE const int32_t& Int(int i) const {return (reinterpret_cast<const int32_t*>(&_data))[i];}

  STK_MATH_FORCE_INLINE uint32_t& UInt(int i) {return (reinterpret_cast<uint32_t*>(&_data))[i];}
  STK_MATH_FORCE_INLINE const uint32_t& UInt(int i) const {return (reinterpret_cast<const uint32_t*>(&_data))[i];}

  float _data[simd::nfloats]; // the "_" means you should try not to use this directly
  // it is made public to avoid function call overhead 
  // and/or so the compiler doesn't have to use up one of
  // inlining depths (usually max inlining depth ~5)
};

} // namespace simd
} // namespace stk

