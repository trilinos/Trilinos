// Copyright 2013 Sandia Corporation, Albuquerque, NM.

// IWYU pragma: private, include <stk_simd/Simd.hpp>

#include <array>
#include "./NoSimdSizes.hpp"

namespace stk {
namespace simd {

struct Double {

  STK_MATH_FORCE_INLINE Double() {}

  template <typename T>
  STK_MATH_FORCE_INLINE Double(const T x, typename std::enable_if<std::is_convertible<T,double>::value, void*>::type=0) {

    for (int i=0; i < ndoubles; ++i) _data[i] = double(x);
  }

  STK_MATH_FORCE_INLINE Double(const Double& x) {

    for (int i=0; i < ndoubles; ++i) _data[i] = x._data[i];  
  }

  STK_MATH_FORCE_INLINE Double& operator= (const Double& x) {

    for (int i=0; i < ndoubles; ++i) _data[i] = x._data[i];
    return *this;
  }

  template <typename T>
  STK_MATH_FORCE_INLINE typename std::enable_if<std::is_convertible<T,double>::value, Double&>::type operator= (const T x) {

    for (int i=0; i < ndoubles; ++i) _data[i] = double(x);
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator+= (const Double& a) {

    for (int i=0; i < ndoubles; ++i) _data[i] += a[i];
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator-= (const Double& a) {

    for (int i=0; i < ndoubles; ++i) _data[i] -= a[i];
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator*= (const Double& a) {

    for (int i=0; i < ndoubles; ++i) _data[i] *= a[i];
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator/= (const Double& a) {

    for (int i=0; i < ndoubles; ++i) _data[i] /= a[i];
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator+= (const double a) {

    for (int i=0; i < ndoubles; ++i) _data[i] += a;
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator-= (const double a) {

    for (int i=0; i < ndoubles; ++i) _data[i] -= a;
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator*= (const double a) {

    for (int i=0; i < ndoubles; ++i) _data[i] *= a;
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator/= (const double a) {

    for (int i=0; i < ndoubles; ++i) _data[i] /= a;
    return *this;
  }

  STK_MATH_FORCE_INLINE Double operator-() const {

    Double tmp;
    for (int i=0; i < ndoubles; ++i) tmp[i] = -_data[i];
    return tmp;
  }

  STK_MATH_FORCE_INLINE double& operator[](int i) {return _data[i];}
  STK_MATH_FORCE_INLINE const double& operator[](int i) const {return _data[i];}

  double _data[simd::ndoubles];
};

} // namespace simd
} // namespace stk

