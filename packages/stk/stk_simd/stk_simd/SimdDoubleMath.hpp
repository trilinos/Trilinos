// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// IWYU pragma: private; include "Simd.hpp"

#ifndef STK_INCLUDE_ONLY_STK_SIMD_HEADER
static_assert(false, "Do not include simd impl files directly. Only include stk_simd/Simd.hpp");
#endif

#ifndef STK_SIMD_DOUBLEMATH_HPP
#define STK_SIMD_DOUBLEMATH_HPP

#if (defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)) && !defined(USE_STK_SIMD_NONE)
  #define STK_SIMD_INTEL_ENABLED 1
#else
  #define STK_SIMD_INTEL_ENABLED 0
#endif

namespace stk {
namespace math {

namespace hidden {
  template <typename Func>
  STK_MATH_FORCE_INLINE simd::Double apply_per_lane(const simd::Double& x, Func&& func) {
#if defined(STK_SIMD_USE_NON_UB_METHOD)
    alignas(alignof(simd::Double)) double  in[simd::ndoubles];
    alignas(alignof(simd::Double)) double out[simd::ndoubles];
    simd::store(in, x);
    for (int i = 0; i < simd::ndoubles; ++i) {
      out[i] = func(in[i]);
    }
    return simd::load(out);
#else
    simd::Double tmp;
    for (int i = 0; i < simd::ndoubles; ++i) {
      tmp[i] = func(x[i]);
    }
    return tmp;
#endif
  }

  template <typename Func>
  STK_MATH_FORCE_INLINE simd::Double apply_per_lane(const simd::Double& x, const simd::Double& y, Func&& func) {
#if defined(STK_SIMD_USE_NON_UB_METHOD)
    alignas(alignof(simd::Double)) double in1[simd::ndoubles];
    alignas(alignof(simd::Double)) double in2[simd::ndoubles];
    alignas(alignof(simd::Double)) double out[simd::ndoubles];

    simd::store(in1, x);
    simd::store(in2, y);

    for (int i = 0; i < simd::ndoubles; ++i) {
      out[i] = func(in1[i], in2[i]);
    }
    return simd::load(out);
#else
    simd::Double tmp;
    for (int i = 0; i < simd::ndoubles; ++i) {
      tmp[i] = func(x[i], y[i]);
    }
    return tmp;
#endif
  }

} // end hidden namespace

STK_MATH_FORCE_INLINE simd::Double fmadd(const simd::Double& a, const simd::Double& b, const simd::Double& c) {
  return simd::Double(fma(a._data, b._data, c._data));
}

STK_MATH_FORCE_INLINE simd::Double sqrt(const simd::Double& x) {
  return simd::Double(SIMD_NAMESPACE::sqrt(x._data));
}

STK_MATH_FORCE_INLINE simd::Double cbrt(const simd::Double& x) {
#if STK_SIMD_INTEL_ENABLED
  return simd::Double(SIMD_NAMESPACE::cbrt(x._data));
#else
  return hidden::apply_per_lane(x, [](double val) { return std::cbrt(val); });
#endif
}

STK_MATH_FORCE_INLINE simd::Double log(const simd::Double& x) {
#if STK_SIMD_INTEL_ENABLED
  return simd::Double(SIMD_NAMESPACE::log(x._data));
#else
  return hidden::apply_per_lane(x, [](double val) { return std::log(val); });
#endif
}

STK_MATH_FORCE_INLINE simd::Double log10(const simd::Double& x) {
  return hidden::apply_per_lane(x, [](double val) { return std::log10(val); });
}

STK_MATH_FORCE_INLINE simd::Double exp(const simd::Double& x) {
#if STK_SIMD_INTEL_ENABLED
  return simd::Double(SIMD_NAMESPACE::exp(x._data));
#else
  return hidden::apply_per_lane(x, [](double val) { return std::exp(val); });
#endif
}



/**
 * @brief Computes the power of a SIMD double raised to an integer exponent.
 *
 * Uses exponentiation by squaring to efficiently compute x^y. If the 
 * exponent y is negative, the result is computed as the reciprocal of x^|y|.
 */
STK_MATH_FORCE_INLINE simd::Double pow(const simd::Double& x, const int y) {
  simd::Double base = x;
  simd::Double result = simd::Double(1.0);
  int n = y;
  const bool negativeExponent = (n < 0);
  if (negativeExponent) n = -n;
  while (n) {
    if (n & 1) result *= base;
    base *= base;
    n >>= 1;
  }
  return negativeExponent ? simd::Double(1.0) / result : result;
}

STK_MATH_FORCE_INLINE simd::Double sin(const simd::Double& a) {
  return hidden::apply_per_lane(a, [](double val) { return std::sin(val); });
}

STK_MATH_FORCE_INLINE simd::Double cos(const simd::Double& a) {
  return hidden::apply_per_lane(a, [](double val) { return std::cos(val); });
}

STK_MATH_FORCE_INLINE simd::Double tan(const simd::Double& a) {
  return hidden::apply_per_lane(a, [](double val) { return std::tan(val); });
}

STK_MATH_FORCE_INLINE simd::Double sinh(const simd::Double& a) {
  return hidden::apply_per_lane(a, [](double val) { return std::sinh(val); });
}

STK_MATH_FORCE_INLINE simd::Double cosh(const simd::Double& a) {
  return hidden::apply_per_lane(a, [](double val) { return std::cosh(val); });
}

STK_MATH_FORCE_INLINE simd::Double tanh(const simd::Double& a) {
  return hidden::apply_per_lane(a, [](double val) { return std::tanh(val); });
}

// Inverse trigonometric functions
STK_MATH_FORCE_INLINE simd::Double asin(const simd::Double& a) {
  return hidden::apply_per_lane(a, [](double val) { return std::asin(val); });
}

STK_MATH_FORCE_INLINE simd::Double acos(const simd::Double& a) {
  return hidden::apply_per_lane(a, [](double val) { return std::acos(val); });
}

STK_MATH_FORCE_INLINE simd::Double atan(const simd::Double& a) {
  return hidden::apply_per_lane(a, [](double val) { return std::atan(val); });
}

STK_MATH_FORCE_INLINE simd::Double atan2(const simd::Double& a, const simd::Double& b) {
  return hidden::apply_per_lane(a, b, [](double A, double B) { return std::atan2(A, B); });
}

STK_MATH_FORCE_INLINE simd::Double asinh(const simd::Double& a) {
  return hidden::apply_per_lane(a, [](double val) { return std::asinh(val); });
}

STK_MATH_FORCE_INLINE simd::Double acosh(const simd::Double& a) {
  return hidden::apply_per_lane(a, [](double val) { return std::acosh(val); });
}

STK_MATH_FORCE_INLINE simd::Double atanh(const simd::Double& a) {
  return hidden::apply_per_lane(a, [](double val) { return std::atanh(val); });
}

STK_MATH_FORCE_INLINE simd::Double erf(const simd::Double& a) {
  return hidden::apply_per_lane(a, [](double val) { return std::erf(val); });
}

STK_MATH_FORCE_INLINE simd::Double multiplysign(const simd::Double& x, const simd::Double& y) {
  return simd::Double(SIMD_NAMESPACE::multiplysign(x._data, y._data));
}

STK_MATH_FORCE_INLINE simd::Double copysign(const simd::Double& x, const simd::Double& y) {
  return simd::Double(SIMD_NAMESPACE::copysign(x._data, y._data));
}

STK_MATH_FORCE_INLINE simd::Double abs(const simd::Double& x) {
  return simd::Double(SIMD_NAMESPACE::abs(x._data));
}

STK_MATH_FORCE_INLINE simd::Double min(const simd::Double& x, const simd::Double& y) {
  return simd::Double(SIMD_NAMESPACE::min(x._data, y._data));
}

STK_MATH_FORCE_INLINE simd::Double max(const simd::Double& x, const simd::Double& y) {
  return simd::Double(SIMD_NAMESPACE::max(x._data, y._data));
}

STK_MATH_FORCE_INLINE simd::Bool isnan(const simd::Double& a) {
  // Relies on IEEE property: NaN != NaN
  return a != a;
}

STK_MATH_FORCE_INLINE simd::Double if_then_else(const simd::Bool& b, const simd::Double& v1, const simd::Double& v2) {
  return simd::Double(SIMD_NAMESPACE::choose(b._data, v1._data, v2._data));
}

STK_MATH_FORCE_INLINE simd::Double if_then_else_zero(const simd::Bool& b, const simd::Double& v) {
  return simd::Double(SIMD_NAMESPACE::choose(b._data, v._data, SIMD_NAMESPACE::native_simd<double>(0.0)));
}

template <typename T>
STK_MATH_FORCE_INLINE bool are_all_equal(const T& vec, double value) {
  for (int i = 0; i < stk::simd::ndoubles; ++i) {
    if (vec[i] != value) {
      return false;
    }
  }
  return true;
}

#if STK_SIMD_INTEL_ENABLED
/**
 * @brief Compute SIMD element-wise pow(x, y) = x^y using vectorized log/exp.
 *
 * This is the "basic" fast version: it assumes well-behaved inputs and does not
 * fully replicate all IEEE corner cases (like std::pow). It is intended for
 * performance-sensitive kernels where speed matters more than covering every
 * special case.
 *
 * Core idea:
 *   x^y = exp(y * log(x))     when x > 0
 *
 * When x < 0 and y is an integer:
 *   - The magnitude is computed as |x|^y
 *   - The sign is determined by whether y is odd (result negative) or even
 *     (result positive)
 *
 * Vectorization:
 *   - Uses Intel SVML (Short Vector Math Library) intrinsics for exp/log,
 *     which compute across all SIMD lanes in parallel.
 *   - This avoids looping over each element and gives ~2.5x speedup
 *     compared to scalar pow() per lane.
 *
 * Limitations:
 *   - If x < 0 and y is not an integer, this returns a value instead of NaN.
 *   - No explicit handling of x == 0, y == 0, infinities, or NaN inputs.
 *
 * References:
 *   - "Elementary Functions: Algorithms and Implementation" by Jean-Michel Muller
 *   - SVML intrinsics documentation (Intel oneAPI Math Kernel Library)
 *   - IEEE-754 standard for floating-point arithmetic
 */
STK_MATH_FORCE_INLINE simd::Double pow_vectorized(const simd::Double& x, const simd::Double& y) {
  using namespace stk::math;

  // Quick return for special cases
  if (are_all_equal(y, 0.0)) {
    return simd::Double(1.0); // x^0 = 1 for all lanes
  }
  if (are_all_equal(y, 1.0)) {
    return x; // x^1 = x for all lanes
  }
  if (are_all_equal(y, -1.0)) {
    return simd::Double(1.0) / x; // x^-1 = 1 / x for all lanes
  }

  // Step 1. Compute the magnitude
  //
  // For any real x > 0, we can compute x^y as:
  //    magnitude = exp(y * log(x))
  //
  // If x < 0, we use |x| in the magnitude, and later adjust the sign.
  const simd::Double ax  = abs(x);      // |x|
  const simd::Double t   = y * log(ax); // y * log(|x|)
  const simd::Double mag = exp(t);      // exp(y * log(|x|)) = |x|^y

  // Step 2. Determine if exponent y is odd or even (integer case only)
  //
  // Trick: Divide y by 2, round to nearest integer, and check the fractional part.
  // If (y/2 - round(y/2)) == 0.5, then y is odd.
  //
  // We use a "round to nearest" trick for doubles: add a large constant (2^52+2^51),
  // then subtract it back. This works because doubles have 52 bits of precision
  // in the mantissa. See: "Hacker's Delight" (Henry S. Warren, 2002).
  auto round_nearest_abs = [](const simd::Double& v_abs) {
    constexpr double magic = 6755399441055744.0; // 2^52 + 2^51
    return (v_abs + simd::Double(magic)) - simd::Double(magic);
  };

  const simd::Double one(1.0);
  const simd::Double half(0.5);
  const simd::Double yh_abs   = abs(y) * half; // |y| / 2
  const simd::Double yh_rnabs = round_nearest_abs(yh_abs);
  const auto         odd      = (abs(yh_abs - yh_rnabs) == half);

  // If y is odd, the sign is -1; otherwise, the sign is +1.
  const simd::Double sign_if_neg = if_then_else(odd, simd::Double(-1.0), one);

  // 3) Form the lane-wise sign and return with a single multiply
  // sign = +1 when x >= 0; otherwise ±1 depending on parity.
  const simd::Double sign = if_then_else(x < simd::Double(0.0), sign_if_neg, one);

  return sign * mag;
}
#endif

STK_MATH_FORCE_INLINE simd::Double pow(const simd::Double& x, const simd::Double& y) {
#if STK_SIMD_INTEL_ENABLED
  return pow_vectorized(x, y);
#else
  return hidden::apply_per_lane(x, y, [](double X, double Y) { return std::pow(X, Y); });
#endif
}

STK_MATH_FORCE_INLINE simd::Double pow(const simd::Double& x, const double y) {
#if STK_SIMD_INTEL_ENABLED
  return pow(x, simd::Double(y));
#else
  return hidden::apply_per_lane(x, [y](double base) { return std::pow(base, y); });
#endif
}

} // namespace math
} // namespace stk

#endif //STK_SIMD_DOUBLEMATH_HPP

