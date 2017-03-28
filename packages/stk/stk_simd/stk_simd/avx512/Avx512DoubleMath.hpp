
namespace stk {
namespace math {

namespace hidden {
inline static double CBRT(const double x) { return cbrt(x); }
static const simd::Double SIGN_MASK(-0.0);
static const simd::Double ZERO(0.0);
}

inline simd::Double fmadd(const simd::Double& a, const simd::Double& b, const simd::Double& c) {
  return simd::Double(_mm512_fmadd_pd(a._data, b._data, c._data));
}

inline simd::Double sqrt(const simd::Double& x) {
  return simd::Double(_mm512_sqrt_pd(x._data));
}
  
inline simd::Double cbrt(const simd::Double& x) {
#if defined(__INTEL_COMPILER)
  return simd::Double(_mm512_cbrt_pd(x._data));
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
  return simd::Double(_mm512_log_pd(x._data));
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
  return simd::Double(_mm512_exp_pd(x._data));
#else
  simd::Double tmp;
  for (int n=0; n < simd::ndoubles; ++n) {
    tmp[n] = std::exp(x[n]);
  }
  return tmp;
#endif
}

inline simd::Double pow(const simd::Double& x, const simd::Double& y) {
#if defined(__INTEL_COMPILER)
  return simd::Double(_mm512_pow_pd(x._data, y._data));
#else
  simd::Double tmp;
  for (int n=0; n < simd::ndoubles; ++n) {
    tmp[n] = std::pow(x[n],y[n]);
  }
  return tmp;
#endif
}

inline simd::Double pow(const simd::Double& x, double y) {
#if defined(__INTEL_COMPILER)
  return simd::Double(_mm512_pow_pd(x._data, _mm512_set1_pd(y)));
#else
  simd::Double tmp;
  for (int n=0; n < simd::ndoubles; ++n) {
    tmp[n] = std::pow(x[n],y);
  }
  return tmp;
#endif
}

inline simd::Double pow(const simd::Double& x, int y) {
  simd::Double tmp;
  for (int n=0; n < simd::ndoubles; ++n) {
    tmp[n] = std::pow(x[n],y);
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
  __m512i sign_y =  _mm512_and_epi64(reinterpret_cast<const __m512i&>(hidden::SIGN_MASK._data),reinterpret_cast<const __m512i&>(y._data));
  __m512i sol = _mm512_xor_epi64(reinterpret_cast<const __m512i&>(x._data),sign_y);
  return simd::Double( reinterpret_cast<__m512d&>(sol) );
}
  
inline simd::Double copysign(const simd::Double& x, const simd::Double& y) { // return abs(x) times sign of y
  __m512i pos_x = _mm512_andnot_epi64(reinterpret_cast<const __m512i&>(hidden::SIGN_MASK._data), reinterpret_cast<const __m512i&>(x._data));
  __m512i sign_y = _mm512_and_epi64(reinterpret_cast<const __m512i&>(hidden::SIGN_MASK._data),reinterpret_cast<const __m512i&>(y._data));
  __m512i sol = _mm512_xor_epi64(pos_x,sign_y);
  return simd::Double( reinterpret_cast<__m512d&>(sol) );
}

inline simd::Double abs(const simd::Double& x) {
  __m512i sol = _mm512_andnot_epi64(  reinterpret_cast<const __m512i&>(hidden::SIGN_MASK._data), reinterpret_cast<const __m512i&>(x._data) ); 
  return simd::Double( reinterpret_cast<__m512d&>(sol) );
}

inline simd::Double min(const simd::Double& x, const simd::Double& y) {
  return simd::Double(_mm512_min_pd(x._data, y._data));
}

inline simd::Double max(const simd::Double& x, const simd::Double& y) {
  return simd::Double(_mm512_max_pd(x._data, y._data));
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
  return simd::Double( _mm512_mask_add_pd(v2._data, b._data, v1._data, hidden::ZERO._data) );
}

inline simd::Double if_then_else_zero(const simd::Bool& b, const simd::Double& v) {
  return simd::Double( _mm512_mask_add_pd(hidden::ZERO._data, b._data, v._data, hidden::ZERO._data) );
}

inline simd::Double if_not_then_else(const simd::Bool& b, const simd::Double& v1, const simd::Double& v2) {
  return simd::Double( _mm512_mask_add_pd(v1._data, b._data, v2._data, hidden::ZERO._data) );
}

inline simd::Double if_not_then_else_zero(const simd::Bool& b, const simd::Double& v) {
  return simd::Double( _mm512_mask_add_pd(v._data, b._data, hidden::ZERO._data, hidden::ZERO._data) );
}

} // namespace math
} // namespace stk
