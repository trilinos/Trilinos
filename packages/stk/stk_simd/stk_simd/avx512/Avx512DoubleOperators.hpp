
namespace stk {
namespace simd {

namespace hidden {
static const simd::Bool TRUE_VAL(true);
}

inline simd::Double operator+ (const simd::Double& a, const simd::Double& b) {
  return simd::Double(_mm512_add_pd(a._data, b._data));
}

inline simd::Double operator- (const simd::Double& a, const simd::Double& b) {
  return simd::Double(_mm512_sub_pd(a._data, b._data));
}

inline simd::Double operator* (const simd::Double& a, const simd::Double& b) {
  return simd::Double(_mm512_mul_pd(a._data, b._data));
}
   
inline simd::Double operator/ (const simd::Double& a, const simd::Double& b) {
  return simd::Double(_mm512_div_pd(a._data, b._data));
}

inline simd::Double operator+ (const double a, const simd::Double& b) {
  return simd::Double(_mm512_add_pd(_mm512_set1_pd(a),b._data));
}

inline simd::Double operator+ (const simd::Double& a, const double b) {
  return simd::Double(_mm512_add_pd(a._data, _mm512_set1_pd(b)));
}

inline simd::Double operator- (const double a, const simd::Double& b) {
  return simd::Double(_mm512_sub_pd(_mm512_set1_pd(a),b._data));
}
  
inline simd::Double operator- (const simd::Double& a, const double b) {
  return simd::Double(_mm512_sub_pd(a._data, _mm512_set1_pd(b)));
}

inline simd::Double operator* (const double a, const simd::Double& b) {
  return simd::Double(_mm512_mul_pd(_mm512_set1_pd(a),b._data));
}

inline simd::Double operator* (const simd::Double& a, const double b) {
  return simd::Double(_mm512_mul_pd(a._data, _mm512_set1_pd(b)));
}

inline simd::Double operator/ (const double a, const simd::Double& b) {
  return simd::Double(_mm512_div_pd(_mm512_set1_pd(a),b._data));
}
  
inline simd::Double operator/ (const simd::Double& a, const double b) {
  return simd::Double(_mm512_div_pd(a._data, _mm512_set1_pd(b)));
}

inline simd::Bool operator< (const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm512_cmplt_pd_mask(x._data, y._data));
}
    
inline simd::Bool operator<=(const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm512_cmple_pd_mask(x._data, y._data));
}

inline simd::Bool operator> (const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm512_cmplt_pd_mask(y._data, x._data));
}

inline simd::Bool operator>=(const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm512_cmple_pd_mask(y._data, x._data));
}

inline simd::Bool operator==(const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm512_cmpeq_pd_mask(x._data, y._data));
}
  
inline simd::Bool operator!=(const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm512_cmpneq_pd_mask(x._data, y._data));
}

inline simd::Bool operator&& (const simd::Bool& x, const simd::Bool& y) {
  return simd::Bool( (__mmask8)_mm512_kand(x._data, y._data));
}

inline simd::Bool operator|| (const simd::Bool& x, const simd::Bool& y) {
  return simd::Bool( (__mmask8)_mm512_kor(x._data, y._data));
}

inline simd::Bool operator! (const simd::Bool& x) {
  return simd::Bool( (__mmask8)_mm512_kxor(hidden::TRUE_VAL._data, x._data) );
}

} // namespace simd
} // namespace stk
