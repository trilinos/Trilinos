
namespace stk {
namespace simd {

namespace hidden {
static const simd::Boolf TRUE_VALf(true);
}

STK_MATH_FORCE_INLINE simd::Float operator+ (const simd::Float& a, const simd::Float& b) {
  return simd::Float(_mm256_add_ps(a._data, b._data));
}

STK_MATH_FORCE_INLINE simd::Float operator- (const simd::Float& a, const simd::Float& b) {
  return simd::Float(_mm256_sub_ps(a._data, b._data));
}

STK_MATH_FORCE_INLINE simd::Float operator* (const simd::Float& a, const simd::Float& b) {
  return simd::Float(_mm256_mul_ps(a._data, b._data));
}
   
STK_MATH_FORCE_INLINE simd::Float operator/ (const simd::Float& a, const simd::Float& b) {
  return simd::Float(_mm256_div_ps(a._data, b._data));
}

STK_MATH_FORCE_INLINE simd::Float operator+ (const float a, const simd::Float& b) {
  return simd::Float(_mm256_add_ps(_mm256_set1_ps(a),b._data));
}

STK_MATH_FORCE_INLINE simd::Float operator+ (const simd::Float& a, const float b) {
  return simd::Float(_mm256_add_ps(a._data, _mm256_set1_ps(b)));
}

STK_MATH_FORCE_INLINE simd::Float operator- (const float a, const simd::Float& b) {
  return simd::Float(_mm256_sub_ps(_mm256_set1_ps(a),b._data));
}
  
STK_MATH_FORCE_INLINE simd::Float operator- (const simd::Float& a, const float b) {
  return simd::Float(_mm256_sub_ps(a._data, _mm256_set1_ps(b)));
}

STK_MATH_FORCE_INLINE simd::Float operator* (const float a, const simd::Float& b) {
  return simd::Float(_mm256_mul_ps(_mm256_set1_ps(a),b._data));
}

STK_MATH_FORCE_INLINE simd::Float operator* (const simd::Float& a, const float b) {
  return simd::Float(_mm256_mul_ps(a._data, _mm256_set1_ps(b)));
}

STK_MATH_FORCE_INLINE simd::Float operator/ (const float a, const simd::Float& b) {
  return simd::Float(_mm256_div_ps(_mm256_set1_ps(a),b._data));
}
  
STK_MATH_FORCE_INLINE simd::Float operator/ (const simd::Float& a, const float b) {
  return simd::Float(_mm256_div_ps(a._data, _mm256_set1_ps(b)));
}

STK_MATH_FORCE_INLINE simd::Boolf operator< (const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm256_cmp_ps(x._data, y._data, _CMP_LT_OQ));
}

STK_MATH_FORCE_INLINE simd::Boolf operator<=(const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm256_cmp_ps(x._data, y._data, _CMP_LE_OQ));
}

STK_MATH_FORCE_INLINE simd::Boolf operator> (const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm256_cmp_ps(y._data, x._data, _CMP_LT_OQ));
}

STK_MATH_FORCE_INLINE simd::Boolf operator>=(const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm256_cmp_ps(y._data, x._data, _CMP_LE_OQ));
}

STK_MATH_FORCE_INLINE simd::Boolf operator==(const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm256_cmp_ps(x._data, y._data, _CMP_EQ_OQ));
}
  
STK_MATH_FORCE_INLINE simd::Boolf operator!=(const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm256_cmp_ps(x._data, y._data, _CMP_NEQ_UQ));
}

STK_MATH_FORCE_INLINE simd::Boolf operator&& (const simd::Boolf& x, const simd::Boolf& y) {
  return simd::Boolf(_mm256_and_ps(x._data, y._data));
}

STK_MATH_FORCE_INLINE simd::Boolf operator|| (const simd::Boolf& x, const simd::Boolf& y) {
  return simd::Boolf(_mm256_or_ps(x._data, y._data));
}

STK_MATH_FORCE_INLINE simd::Boolf operator! (const simd::Boolf& x) {
  return simd::Boolf( _mm256_xor_ps(simd::TRUE_VALf._data, x._data) );
}

} // namespace simd
} // namespace stk
