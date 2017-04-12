
#ifndef STK_MATH_H
#define STK_MATH_H

#include <stdio.h>
#include <cmath>
#include <assert.h>
#include <Kokkos_Macros.hpp>

#define STK_MATH_FORCE_INLINE KOKKOS_FORCEINLINE_FUNCTION

namespace stk {
namespace math {

namespace hidden {
static STK_MATH_FORCE_INLINE double Cbrt(const double x) { return cbrt(x); }
static STK_MATH_FORCE_INLINE float CbrtF(const float x) { return cbrtf(x); }
}

// double

STK_MATH_FORCE_INLINE double fmadd(const double a, const double b, const double c) {
  return a*b+c;
}

STK_MATH_FORCE_INLINE double sqrt(const double x) {
  return std::sqrt(x);
}

STK_MATH_FORCE_INLINE double cbrt(const double x) {
  return hidden::Cbrt(x);
}

STK_MATH_FORCE_INLINE double log(const double x) {
  return std::log(x);
}

STK_MATH_FORCE_INLINE double exp(const double x) {
  return std::exp(x);
}
   
STK_MATH_FORCE_INLINE double pow(const double x, const double y) {
  return std::pow(x, y);
}

STK_MATH_FORCE_INLINE double pow(const double x, const int y) {
  return std::pow(x, y);
}

STK_MATH_FORCE_INLINE double sin(const double a) {
  return std::sin(a);
}

STK_MATH_FORCE_INLINE double cos(const double a) {
  return std::cos(a);
}

STK_MATH_FORCE_INLINE double tan(const double a) {
  return std::tan(a);
}

STK_MATH_FORCE_INLINE double asin(const double a) {
  return std::asin(a);
}

STK_MATH_FORCE_INLINE double acos(const double a) {
  return std::acos(a);
}

STK_MATH_FORCE_INLINE double atan(const double a) {
  return std::atan(a);
}

STK_MATH_FORCE_INLINE double atan2(const double a, const double b) {
  return std::atan2(a, b);
}

STK_MATH_FORCE_INLINE double multiplysign(const double x, const double y) { // return x times sign of y
  return x * ( 1 - 2*(y<0) );
}

STK_MATH_FORCE_INLINE double copysign(const double x, const double y) { // return abs(x) times sign of y
  return std::abs(x) * ( 1 - 2*(y<0) );
}

STK_MATH_FORCE_INLINE double abs(const double x) {
  return std::abs(x);
}

STK_MATH_FORCE_INLINE double max(const double x, const double y) {
  return x > y ? x : y;
}

STK_MATH_FORCE_INLINE double min(const double x, const double y) {
  return x < y ? x : y;
}

STK_MATH_FORCE_INLINE bool isnan(const double a) {
  return a!=a;
}

STK_MATH_FORCE_INLINE double if_then_else(const bool b, const double v1, const double v2) {
  return b ? v1 : v2;
}

STK_MATH_FORCE_INLINE double if_then_else_zero(const bool b, const double v) {
  return b ? v : 0.0;
}

STK_MATH_FORCE_INLINE double if_not_then_else(const bool b, const double v1, const double v2) {
  return b ? v2 : v1;
}

STK_MATH_FORCE_INLINE double if_not_then_else_zero(const bool b, const double v) {
  return b ? 0.0 : v;
}

// float

STK_MATH_FORCE_INLINE float fmadd(const float a, const float b, const float c) {
  return a*b+c;
}

STK_MATH_FORCE_INLINE float sqrt(const float x) {
  return std::sqrt(x);
}

STK_MATH_FORCE_INLINE float cbrt(const float x) {
  return hidden::CbrtF(x);
}

STK_MATH_FORCE_INLINE float log(const float x) {
  return std::log(x);
}

STK_MATH_FORCE_INLINE float exp(const float x) {
  return std::exp(x);
}
   
STK_MATH_FORCE_INLINE float pow(const float x, const float y) {
  return std::pow(x, y);
}

STK_MATH_FORCE_INLINE float pow(const float x, const int y) {
  return std::pow(x, y);
}

STK_MATH_FORCE_INLINE float sin(const float a) {
  return std::sin(a);
}

STK_MATH_FORCE_INLINE float cos(const float a) {
  return std::cos(a);
}

STK_MATH_FORCE_INLINE float tan(const float a) {
  return std::tan(a);
}

STK_MATH_FORCE_INLINE float asin(const float a) {
  return std::asin(a);
}

STK_MATH_FORCE_INLINE float acos(const float a) {
  return std::acos(a);
}

STK_MATH_FORCE_INLINE float atan(const float a) {
  return std::atan(a);
}

STK_MATH_FORCE_INLINE float atan2(const float a, const float b) {
  return std::atan2(a,b);
}

STK_MATH_FORCE_INLINE float multiplysign(const float x, const float y) { // return x times sign of y
  return x * ( 1 - 2*(y<0) );
}

STK_MATH_FORCE_INLINE float copysign(const float x, const float y) { // return abs(x) times sign of y
  return std::abs(x) * ( 1 - 2*(y<0) );
}

STK_MATH_FORCE_INLINE float abs(const float x) {
  return std::abs(x);
}

STK_MATH_FORCE_INLINE float max(const float x, const float y) {
  return x > y ? x : y;
}

STK_MATH_FORCE_INLINE float min(const float x, const float y) {
  return x < y ? x : y;
}

STK_MATH_FORCE_INLINE bool isnan(const float a) {
  return a!=a;
}

STK_MATH_FORCE_INLINE float if_then_else(const bool b, const float v1, const float v2) {
  return b ? v1 : v2;
}

STK_MATH_FORCE_INLINE float if_then_else_zero(const bool b, const float v) {
  return b ? v : 0.0;
}

STK_MATH_FORCE_INLINE float if_not_then_else(const bool b, const float v1, const float v2) {
  return b ? v2 : v1;
}

STK_MATH_FORCE_INLINE float if_not_then_else_zero(const bool b, const float v) {
  return b ? 0.0 : v;
}

} // namespace math
} // namespace stk

#endif
