
namespace stk {
namespace math {

namespace hidden {
inline static double CBRT(const double x) { return cbrt(x); }
static const simd::Double SIGN_MASK(-0.0);
}

inline simd::Double fmadd(const simd::Double& a, const simd::Double& b, const simd::Double& c) {
  return simd::Double(_mm_add_pd(_mm_mul_pd(a._data, b._data), c._data));
}

inline simd::Double sqrt(const simd::Double& x) {
  return simd::Double(_mm_sqrt_pd(x._data));
}
  
inline simd::Double cbrt(const simd::Double& x) {
#if defined(__INTEL_COMPILER)
  return simd::Double(_mm_cbrt_pd(x._data));
#else
  simd::Double tmp;
  for (int n=0; n < simd::ndoubles; ++n) {
    tmp[n] = hidden::CBRT(x[n]);
  }
  return tmp;
#endif
}

inline simd::Double log(const simd::Double& x) {
#if defined(__INTEL_COMPILER)
  return simd::Double(_mm_log_pd(x._data));
#else
  simd::Double tmp;
  for (int n=0; n < simd::ndoubles; ++n) {
    tmp[n] = std::log(x[n]);
  }
  return tmp;
#endif
}

inline simd::Double exp(const simd::Double& x) {
#if defined(__INTEL_COMPILER)
  return simd::Double(_mm_exp_pd(x._data));
#else
  simd::Double tmp;
  for (int n=0; n < simd::ndoubles; ++n) {
    tmp[n] = std::exp(x[n]);
  }
  return tmp;
#endif
}

inline simd::Double pow(const simd::Double& x, const simd::Double& y) {
  simd::Double tmp;
  for (int n=0; n < simd::ndoubles; ++n) {
    tmp[n] = std::pow(x[n], y[n]);
  }
  return tmp;
}

inline simd::Double pow(const simd::Double& x, const double y) {
  simd::Double tmp;
  for (int n=0; n < simd::ndoubles; ++n) {
    tmp[n] = std::pow(x[n], y);
  }
  return tmp;
}

inline simd::Double pow(const simd::Double& x, const int y) {
  simd::Double tmp;
  for (int n=0; n < simd::ndoubles; ++n) {
    tmp[n] = std::pow(x[n], y);
  }
  return tmp;
}

inline simd::Double sin(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) {
    tmp[i] = std::sin(a[i]);
  }
  return tmp;
}
 
inline simd::Double cos(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) {
    tmp[i] = std::cos(a[i]);
  }
  return tmp;
}

inline simd::Double tan(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) {
    tmp[i] = std::tan(a[i]);
  }
  return tmp;
}

inline simd::Double asin(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) {
    tmp[i] = std::asin(a[i]);
  }
  return tmp;
}

inline simd::Double acos(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) {
    tmp[i] = std::acos(a[i]);
  }
  return tmp;
}

inline simd::Double atan(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) {
    tmp[i] = std::atan(a[i]);
  }
  return tmp;
}

inline simd::Double atan2(const simd::Double& a, const simd::Double& b) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) {
    tmp[i] = std::atan2(a[i],b[i]);
  }
  return tmp;
}

inline simd::Double multiplysign(const simd::Double& x, const simd::Double& y) { // return x times sign of y
  return simd::Double(_mm_xor_pd(x._data, _mm_and_pd(hidden::SIGN_MASK._data, y._data)));
}
  
inline simd::Double copysign(const simd::Double& x, const simd::Double& y) { // return abs(x) times sign of y
  return simd::Double(_mm_xor_pd(_mm_andnot_pd(hidden::SIGN_MASK._data, x._data), _mm_and_pd(hidden::SIGN_MASK._data,y._data)));
}

inline simd::Double abs(const simd::Double& x) {
  return simd::Double(_mm_andnot_pd(hidden::SIGN_MASK._data, x._data)); // !sign_mask & x
}

inline simd::Double min(const simd::Double& x, const simd::Double& y) {
  return simd::Double(_mm_min_pd(x._data, y._data));
}

inline simd::Double max(const simd::Double& x, const simd::Double& y) {
  return simd::Double(_mm_max_pd(x._data, y._data));
}

inline simd::Bool isnan(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) {
    if ( a[i] != a[i] ) {
      tmp[i] = 1.0;
    } else {
      tmp[i] = 0.0;
    }
  }
  return tmp > simd::Double(0.5f);
}

inline simd::Double if_then_else(const simd::Bool& b, const simd::Double& v1, const simd::Double& v2) {
  return simd::Double( _mm_add_pd(_mm_and_pd(b._data, v1._data), _mm_andnot_pd(b._data, v2._data)) );
}

inline simd::Double if_then_else_zero(const simd::Bool& b, const simd::Double& v) {
  return simd::Double( _mm_and_pd(b._data, v._data) );
}

inline simd::Double if_not_then_else(const simd::Bool& b, const simd::Double& v1, const simd::Double& v2) {
  return simd::Double( _mm_add_pd(_mm_and_pd(b._data, v2._data), _mm_andnot_pd(b._data, v1._data)) );
}

inline simd::Double if_not_then_else_zero(const simd::Bool& b, const simd::Double& v) {
  return simd::Double( _mm_andnot_pd(b._data, v._data) );
}

} // namespace math
} // namespace stk
