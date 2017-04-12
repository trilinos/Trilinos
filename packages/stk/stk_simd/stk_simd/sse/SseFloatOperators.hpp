
namespace stk {
namespace simd {

namespace hidden {
static const simd::Boolf TRUE_VALf(true);
}

inline simd::Float operator+ (const simd::Float& a, const simd::Float& b) {
  return simd::Float(_mm_add_ps(a._data,b._data));
}

inline simd::Float operator- (const simd::Float& a, const simd::Float& b) {
  return simd::Float(_mm_sub_ps(a._data,b._data));
}

inline simd::Float operator* (const simd::Float& a, const simd::Float& b) {
  return simd::Float(_mm_mul_ps(a._data,b._data));
}
   
inline simd::Float operator/ (const simd::Float& a, const simd::Float& b) {
  return simd::Float(_mm_div_ps(a._data,b._data));
}

inline simd::Float operator+ (const float a, const simd::Float& b) {
  return simd::Float(_mm_add_ps(_mm_set1_ps(a),b._data));
}

inline simd::Float operator+ (const simd::Float& a, const float b) {
  return simd::Float(_mm_add_ps(a._data,_mm_set1_ps(b)));
}

inline simd::Float operator- (const float a, const simd::Float& b) {
  return simd::Float(_mm_sub_ps(_mm_set1_ps(a),b._data));
}
  
inline simd::Float operator- (const simd::Float& a, const float b) {
  return simd::Float(_mm_sub_ps(a._data,_mm_set1_ps(b)));
}

inline simd::Float operator* (const float a, const simd::Float& b) {
  return simd::Float(_mm_mul_ps(_mm_set1_ps(a),b._data));
}

inline simd::Float operator* (const simd::Float& a, const float b) {
  return simd::Float(_mm_mul_ps(a._data,_mm_set1_ps(b)));
}

inline simd::Float operator/ (const float a, const simd::Float& b) {
  return simd::Float(_mm_div_ps(_mm_set1_ps(a),b._data));
}
  
inline simd::Float operator/ (const simd::Float& a, const float b) {
  return simd::Float(_mm_div_ps(a._data,_mm_set1_ps(b)));
}

inline simd::Boolf operator< (const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm_cmplt_ps(x._data,y._data));
}

inline simd::Boolf operator<=(const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm_cmple_ps(x._data,y._data));
}

inline simd::Boolf operator> (const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm_cmplt_ps(y._data,x._data));
}

inline simd::Boolf operator>=(const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm_cmple_ps(y._data,x._data));
}

inline simd::Boolf operator==(const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm_cmpeq_ps(x._data,y._data));
}
  
inline simd::Boolf operator!=(const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm_cmpneq_ps(x._data,y._data));
}

inline simd::Boolf operator&& (const simd::Boolf& x, const simd::Boolf& y) {
  return simd::Boolf(_mm_and_ps(x._data,y._data));
}

inline simd::Boolf operator|| (const simd::Boolf& x, const simd::Boolf& y) {
  return simd::Boolf(_mm_or_ps(x._data,y._data));
}

inline simd::Boolf operator! (const simd::Boolf& x) {
  return simd::Boolf( _mm_xor_ps(hidden::TRUE_VALf._data,x._data) );
}

} // namespace simd
} // namespace stk
