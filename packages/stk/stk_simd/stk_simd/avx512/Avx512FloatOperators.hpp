
namespace stk {
namespace simd {

namespace hidden {
static const simd::Boolf TRUE_VALf(true);
}

inline simd::Float operator+ (const simd::Float& a, const simd::Float& b) {
  return simd::Float(_mm512_add_ps(a._data,b._data));
}

inline simd::Float operator- (const simd::Float& a, const simd::Float& b) {
  return simd::Float(_mm512_sub_ps(a._data,b._data));
}

inline simd::Float operator* (const simd::Float& a, const simd::Float& b) {
  return simd::Float(_mm512_mul_ps(a._data,b._data));
}
   
inline simd::Float operator/ (const simd::Float& a, const simd::Float& b) {
  return simd::Float(_mm512_div_ps(a._data,b._data));
}

inline simd::Float operator+ (const float a, const simd::Float& b) {
  return simd::Float(_mm512_add_ps(_mm512_set1_ps(a),b._data));
}

inline simd::Float operator+ (const simd::Float& a, const float b) {
  return simd::Float(_mm512_add_ps(a._data,_mm512_set1_ps(b)));
}

inline simd::Float operator- (const float a, const simd::Float& b) {
  return simd::Float(_mm512_sub_ps(_mm512_set1_ps(a),b._data));
}
  
inline simd::Float operator- (const simd::Float& a, const float b) {
  return simd::Float(_mm512_sub_ps(a._data,_mm512_set1_ps(b)));
}

inline simd::Float operator* (const float a, const simd::Float& b) {
  return simd::Float(_mm512_mul_ps(_mm512_set1_ps(a),b._data));
}

inline simd::Float operator* (const simd::Float& a, const float b) {
  return simd::Float(_mm512_mul_ps(a._data,_mm512_set1_ps(b)));
}

inline simd::Float operator/ (const float a, const simd::Float& b) {
  return simd::Float(_mm512_div_ps(_mm512_set1_ps(a),b._data));
}
  
inline simd::Float operator/ (const simd::Float& a, const float b) {
  return simd::Float(_mm512_div_ps(a._data,_mm512_set1_ps(b)));
}

inline simd::Boolf operator< (const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm512_cmplt_ps_mask(x._data,y._data));
}

inline simd::Boolf operator<=(const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm512_cmple_ps_mask(x._data,y._data));
}

inline simd::Boolf operator> (const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm512_cmplt_ps_mask(y._data,x._data));
}

inline simd::Boolf operator>=(const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm512_cmple_ps_mask(y._data,x._data));
}

inline simd::Boolf operator==(const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm512_cmpeq_ps_mask(x._data,y._data));
}
  
inline simd::Boolf operator!=(const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm512_cmpneq_ps_mask(x._data,y._data));
}

inline simd::Boolf operator&& (const simd::Boolf& x, const simd::Boolf& y) {
  return simd::Boolf((__mmask16)_mm512_kand(x._data,y._data));
}

inline simd::Boolf operator|| (const simd::Boolf& x, const simd::Boolf& y) {
  return simd::Boolf((__mmask16)_mm512_kor(x._data,y._data));
}

inline simd::Boolf operator! (const simd::Boolf& x) {
  return simd::Boolf((__mmask16)_mm512_kxor(hidden::TRUE_VALf._data,x._data) );
}

} // namespace simd
} // namespace stk
