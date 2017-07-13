
namespace stk {
namespace simd {

namespace hidden {
static const simd::Bool TRUE_VAL(true);
}

inline simd::Double operator+ (const simd::Double& a, const simd::Double& b) {
  return simd::Double(_mm_add_pd(a._data,b._data));
}

inline simd::Double operator- (const simd::Double& a, const simd::Double& b) {
  return simd::Double(_mm_sub_pd(a._data,b._data));
}

inline simd::Double operator* (const simd::Double& a, const simd::Double& b) {
  return simd::Double(_mm_mul_pd(a._data,b._data));
}
   
inline simd::Double operator/ (const simd::Double& a, const simd::Double& b) {
  return simd::Double(_mm_div_pd(a._data,b._data));
}

inline simd::Double operator+ (const double a, const simd::Double& b) {
  return simd::Double(_mm_add_pd(_mm_set1_pd(a),b._data));
}

inline simd::Double operator+ (const simd::Double& a, const double b) {
  return simd::Double(_mm_add_pd(a._data,_mm_set1_pd(b)));
}

inline simd::Double operator- (const double a, const simd::Double& b) {
  return simd::Double(_mm_sub_pd(_mm_set1_pd(a),b._data));
}
  
inline simd::Double operator- (const simd::Double& a, const double b) {
  return simd::Double(_mm_sub_pd(a._data,_mm_set1_pd(b)));
}

inline simd::Double operator* (const double a, const simd::Double& b) {
  return simd::Double(_mm_mul_pd(_mm_set1_pd(a),b._data));
}

inline simd::Double operator* (const simd::Double& a, const double b) {
  return simd::Double(_mm_mul_pd(a._data,_mm_set1_pd(b)));
}

inline simd::Double operator/ (const double a, const simd::Double& b) {
  return simd::Double(_mm_div_pd(_mm_set1_pd(a),b._data));
}
  
inline simd::Double operator/ (const simd::Double& a, const double b) {
  return simd::Double(_mm_div_pd(a._data,_mm_set1_pd(b)));
}

inline simd::Bool operator< (const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm_cmplt_pd(x._data,y._data));
}

inline simd::Bool operator<=(const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm_cmple_pd(x._data,y._data));
}

inline simd::Bool operator> (const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm_cmplt_pd(y._data,x._data));
}

inline simd::Bool operator>=(const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm_cmple_pd(y._data,x._data));
}

inline simd::Bool operator==(const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm_cmpeq_pd(x._data,y._data));
}
  
inline simd::Bool operator!=(const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm_cmpneq_pd(x._data,y._data));
}

inline simd::Bool operator&& (const simd::Bool& x, const simd::Bool& y) {
  return simd::Bool(_mm_and_pd(x._data,y._data));
}

inline simd::Bool operator|| (const simd::Bool& x, const simd::Bool& y) {
  return simd::Bool(_mm_or_pd(x._data,y._data));
}

inline simd::Bool operator! (const simd::Bool& x) {
  return simd::Bool( _mm_xor_pd(hidden::TRUE_VAL._data,x._data) );
}

} // namespace simd
} // namespace stk
